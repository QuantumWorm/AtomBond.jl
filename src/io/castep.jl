# CASTEP readers: `.cell` input files and `.geom` geometry-optimisation /
# MD trajectories. (The `.cell` writer is in export.jl.)

const BOHR = 0.529177210544          # Å (CODATA 2022)

# Length units CASTEP accepts in a block's optional first line.
const CASTEP_UNITS = Dict("ang" => 1.0, "bohr" => BOHR, "a0" => BOHR, "nm" => 10.0,
                          "cm" => 1e8, "m" => 1e10)

# Species token of a CASTEP line: "Si", "Si:1" / "Si:bulk" (labelled), or an atomic number.
function castep_species(tok::AbstractString)
    s = split(tok, ':')[1]
    z = tryparse(Int, s)
    if z !== nothing
        haskey(ELEMENTS_BY_Z, z) || throw(ArgumentError("unknown atomic number $z"))
        return ELEMENTS_BY_Z[z].symbol
    end
    return normalize_symbol(s)
end

"""
Blocks of a .cell file as Dict(lowercase name => non-empty lines), with
comments (`!`, `#`) removed.
"""
function castep_blocks(io::IO)
    blocks = Dict{String,Vector{String}}()
    current = nothing
    for raw in eachline(io)
        line = strip(split(split(raw, '!')[1], '#')[1])
        isempty(line) && continue
        words = split(line)
        key = lowercase(words[1])
        if key == "%block" && length(words) >= 2
            current = lowercase(words[2])
            blocks[current] = String[]
        elseif key == "%endblock"
            current = nothing
        elseif current !== nothing
            push!(blocks[current], line)
        end
    end
    return blocks
end

# Strip an optional unit line from a block; returns (scale to Å, remaining lines).
function block_unit(lines)
    isempty(lines) && return 1.0, lines
    w = lowercase(split(lines[1])[1])
    haskey(CASTEP_UNITS, w) && length(split(lines[1])) == 1 && return CASTEP_UNITS[w], lines[2:end]
    return 1.0, lines
end

numbers(line, r) = [parse(Float64, x) for x in split(line)[r]]

"""
    read_castep_cell(path or io; title="") -> Crystal

CASTEP `.cell` file: lattice from `LATTICE_CART` or `LATTICE_ABC`, atoms from
`POSITIONS_FRAC` or `POSITIONS_ABS` (units ang, bohr, nm, ...; default Å).
Other blocks and keywords are ignored.
"""
read_castep_cell(path::AbstractString; title=splitext(basename(path))[1]) =
    open(io -> read_castep_cell(io; title), path)

function read_castep_cell(io::IO; title::AbstractString="")
    b = castep_blocks(io)
    L = if haskey(b, "lattice_cart")
        s, lines = block_unit(b["lattice_cart"])
        length(lines) >= 3 || throw(ArgumentError("LATTICE_CART needs three vectors"))
        Lattice((s * Vec3(numbers(lines[k], 1:3)) for k in 1:3)...)
    elseif haskey(b, "lattice_abc")
        s, lines = block_unit(b["lattice_abc"])
        length(lines) >= 2 || throw(ArgumentError("LATTICE_ABC needs lengths and angles"))
        Lattice((s .* numbers(lines[1], 1:3))..., numbers(lines[2], 1:3)...)
    else
        throw(ArgumentError("no LATTICE_CART or LATTICE_ABC block"))
    end
    c = if haskey(b, "positions_frac")
        lines = b["positions_frac"]
        species = [castep_species(split(l)[1]) for l in lines]
        Crystal(L, species, [Vec3(numbers(l, 2:4)) for l in lines]; title)
    elseif haskey(b, "positions_abs")
        s, lines = block_unit(b["positions_abs"])
        species = [castep_species(split(l)[1]) for l in lines]
        crystal_from_cartesian(L, species, [s * Vec3(numbers(l, 2:4)) for l in lines]; title)
    else
        throw(ArgumentError("no POSITIONS_FRAC or POSITIONS_ABS block"))
    end
    # Initial magnetic moments: "Fe 0.0 0.0 0.0 SPIN=2.0" (μB).
    lines = get(b, "positions_frac", get(b, "positions_abs", String[]))
    spins = [(m = match(r"spin\s*[=:]\s*([-+0-9.eE]+)"i, l); m === nothing ? NaN : parse(Float64, m[1]))
             for l in lines if length(split(l)) >= 4 && tryparse(Float64, split(l)[2]) !== nothing]
    length(spins) == natoms(c) && any(!isnan, spins) && (c = with_magmoms(c, spins))
    return castep_constraints(c, b)
end

"""
Apply CELL_CONSTRAINTS and IONIC_CONSTRAINTS. An ionic constraint along a
reciprocal direction fixes that fractional coordinate; an atom with three
independent constraints is fixed completely. Other constraints are ignored.
"""
function castep_constraints(c::Crystal, b)
    cell = c.cell_constraints
    if haskey(b, "cell_constraints") && length(b["cell_constraints"]) >= 2
        cell = Tuple(Int.(round.([numbers(b["cell_constraints"][1], 1:3); numbers(b["cell_constraints"][2], 1:3)])))
    end
    fixed = copy(c.fixed)
    if haskey(b, "ionic_constraints")
        B = inv(c.lattice.matrix)
        vectors = Dict{Int,Vector{Vec3}}()
        for l in b["ionic_constraints"]
            w = split(l)
            length(w) >= 6 || continue
            sp, k = castep_species(w[2]), parse(Int, w[3])
            idx = findall(==(sp), c.species)
            k <= length(idx) || continue
            push!(get!(vectors, idx[k], Vec3[]), Vec3(parse.(Float64, w[4:6])))
        end
        for (i, vs) in vectors
            M = hcat(vs...)
            if rank(M; atol=1e-6) == 3
                fixed[i] = (true, true, true)
                continue
            end
            fixed[i] = ntuple(3) do axis
                n = normalize(Vec3(B[axis, :]))
                fixed[i][axis] || any(v -> abs(abs(dot(normalize(v), n)) - 1) < 1e-6, vs)
            end
        end
    end
    return with_constraints(c; fixed, cell=cell)
end

"""
    read_castep_geom_frames(path or io) -> Vector{Crystal}

All frames of a CASTEP `.geom` (or `.md`) file: cell (`<-- h`) and Cartesian
positions (`<-- R`), converted from Bohr to Å. Each frame's title gives its
step and energy (Hartree in the file, shown in eV).
"""
read_castep_geom_frames(path::AbstractString) = open(io -> read_castep_geom_frames(io; name=basename(path)), path)

function read_castep_geom_frames(io::IO; name::AbstractString="geom")
    frames = Crystal[]
    h, species, pos = Vec3[], Symbol[], Vec3[]
    step, energy = "", NaN
    function finish!()
        isempty(pos) && return
        length(h) == 3 || throw(ArgumentError("frame $(length(frames) + 1): expected 3 cell vectors, got $(length(h))"))
        title = "$name, step $step" * (isnan(energy) ? "" : @sprintf(", E = %.6f eV", energy * 27.211386245981))
        push!(frames, crystal_from_cartesian(Lattice(h...), copy(species), copy(pos); title))
        empty!(h); empty!(species); empty!(pos)
    end
    for raw in eachline(io)
        parts = split(raw, "<--")
        if length(parts) < 2
            w = split(raw)
            # A lone number between frames is the step (or time) of the next frame.
            length(w) == 1 && tryparse(Float64, w[1]) !== nothing && (finish!(); step = w[1])
            continue
        end
        tag = strip(parts[2])
        w = split(parts[1])
        if tag == "E"
            finish!()
            energy = parse(Float64, w[1])
        elseif tag == "h"
            !isempty(pos) && finish!()
            push!(h, BOHR * Vec3(parse.(Float64, w[1:3])))
        elseif tag == "R"
            push!(species, castep_species(w[1]))
            push!(pos, BOHR * Vec3(parse.(Float64, w[3:5])))
        end
    end
    finish!()
    isempty(frames) && throw(ArgumentError("no frames found (expected '<-- h' and '<-- R' lines)"))
    return frames
end

"""
    read_castep_geom(path; frame=:last) -> Crystal

One frame of a CASTEP `.geom` file, by default the last (the relaxed
structure); `frame` may also be an index (1 = starting geometry).
"""
function read_castep_geom(path; frame=:last)
    frames = read_castep_geom_frames(path)
    return frame === :last ? frames[end] : frames[frame]
end

# --- CASTEP output (.castep) ---------------------------------------------------

"""
    read_castep_output(path or io; charges=:mulliken) -> Crystal

CASTEP `.castep` output: the last unit cell ("Real Lattice(A)") and the last
"Cell Contents" fractional coordinates, i.e. the final structure of a
geometry optimization, with atomic charges from the last population
analysis: `charges = :mulliken` (the "Charge (e)" column of "Atomic
Populations (Mulliken)"; Hirshfeld charges when there are no Mulliken ones),
`:hirshfeld`, or `:none`. Magnetic moments come from the Mulliken "Spin
(hbar/2)" column of spin-polarized runs (hbar/2 units = μB).
"""
read_castep_output(path::AbstractString; kw...) =
    open(io -> read_castep_output(io; title=splitext(basename(path))[1], kw...), path)

function read_castep_output(io::IO; title::AbstractString="", charges::Symbol=:mulliken)
    lines = readlines(io)
    # Last lattice: the three lines after "Real Lattice(A)", first three numbers each.
    k = findlast(l -> occursin(r"Real Lattice\s*\(\s*A\s*\)"i, l), lines)
    k === nothing && throw(ArgumentError("no 'Real Lattice(A)' block (not a .castep output?)"))
    L = Lattice((Vec3(parse.(Float64, split(lines[k+i])[1:3])) for i in 1:3)...)
    # Last "Cell Contents" atom table: lines "x  Si   1   0.0 0.0 0.0   x".
    k = findlast(l -> occursin(r"Fractional coordinates of atoms"i, l), lines)
    k === nothing && throw(ArgumentError("no 'Fractional coordinates of atoms' table"))
    species, frac, ion = Symbol[], Vec3[], Dict{Tuple{Symbol,Int},Int}()
    for l in lines[k+1:end]
        w = split(strip(l, ['x', ' ']))
        startswith(strip(l), "xxxx") && !isempty(species) && break
        length(w) >= 5 && tryparse(Int, w[2]) !== nothing && all(t -> tryparse(Float64, t) !== nothing, w[3:5]) || continue
        s = castep_species(w[1])
        push!(species, s)
        push!(frac, Vec3(parse.(Float64, w[3:5])))
        ion[(s, parse(Int, w[2]))] = length(species)
    end
    isempty(species) && throw(ArgumentError("no atoms in the last 'Cell Contents' table"))
    m = fill(NaN, length(species))
    for (s, n, value) in something(castep_population(lines, r"Atomic Populations \(Mulliken\)"i, "Spin"), [])
        i = get(ion, (s, n), 0)
        i > 0 && (m[i] = value)          # spin in hbar/2 = moment in μB (g = 2)
    end
    q = fill(NaN, length(species))
    if charges !== :none
        table = charges === :hirshfeld ? nothing : castep_population(lines, r"Atomic Populations \(Mulliken\)"i, "Charge")
        table === nothing && charges in (:mulliken, :hirshfeld) &&
            (table = castep_population(lines, r"Hirshfeld Analysis"i, "Hirshfeld"))   # "Hirshfeld Charge (e)"
        for (s, n, value) in something(table, [])
            i = get(ion, (s, n), 0)
            i > 0 && (q[i] = value)
        end
    end
    return Crystal(L, species, frac; title, charges=q, magmoms=m)
end

"""
Rows (species, ion, value) of the last population table after a `heading`
line, or `nothing`. The value is the column whose header label starts with
`column`: "Charge" or "Spin" in "Species Ion [Spin] s p d f Total Charge (e)
[Spin (hbar/2)]" (Mulliken), "Hirshfeld" in "Species Ion Hirshfeld Charge (e)".

Layouts differ between CASTEP versions and spin settings (an up/dn column, an
ion split over two lines), so values are matched to the labels from the
right-hand end of the fullest line of each ion.
"""
function castep_population(lines, heading::Regex, column::AbstractString)
    k = findlast(l -> occursin(heading, l), lines)
    k === nothing && return nothing
    h = findnext(l -> occursin(r"^\s*Species\s+Ion"i, l), lines, k)
    h === nothing && return nothing
    # One label per value column: "Hirshfeld Charge" is one label, unit words
    # in brackets are not columns, and a leading "Spin" holds up:/dn:.
    labels = filter(t -> !startswith(t, "("), split(replace(lines[h], r"Hirshfeld\s+Charge" => "Hirshfeld"))[3:end])
    length(labels) > 1 && labels[1] == "Spin" && (labels = labels[2:end])
    col = column == "Spin" ? findlast(t -> startswith(t, "Spin"), labels) : findfirst(t -> startswith(t, column), labels)
    col === nothing && return nothing
    # Rows sit between the "====" lines below the header; an ion may continue
    # on lines without species and number.
    first = findnext(l -> startswith(strip(l), "=="), lines, h)
    first === nothing && return nothing
    groups = Tuple{Symbol,Int,Vector{Vector{Float64}}}[]
    numbers(ws) = [x for x in (tryparse(Float64, t) for t in ws) if x !== nothing]
    for l in lines[first+1:end]
        startswith(strip(l), "==") && break
        w = split(l)
        isempty(w) && continue
        if length(w) >= 2 && tryparse(Int, w[2]) !== nothing && tryparse(Float64, w[1]) === nothing
            push!(groups, (castep_species(w[1]), parse(Int, w[2]), [numbers(w[3:end])]))
        elseif !isempty(groups)
            push!(groups[end][3], numbers(w))
        end
    end
    rows = Tuple{Symbol,Int,Float64}[]
    for (sp, n, candidates) in groups
        nums = argmax(length, candidates)
        i = length(nums) - (length(labels) - col)
        1 <= i <= length(nums) && push!(rows, (sp, n, nums[i]))
    end
    return isempty(rows) ? nothing : rows
end
