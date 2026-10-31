# Inputs of other DFT codes: FHI-aims geometry.in (read and write), SIESTA
# .fdf and ABINIT input (read).

# --- FHI-aims ------------------------------------------------------------------

"""
    read_aims(path or io) -> Crystal

FHI-aims `geometry.in`: `lattice_vector` lines (none for a molecule, which is
boxed), atoms as `atom x y z El` (Å) or `atom_frac f1 f2 f3 El`.
"""
read_aims(path::AbstractString) = open(io -> read_aims(io; title=splitext(basename(path))[1]), path)

function read_aims(io::IO; title::AbstractString="")
    vecs, species, cart, isfrac = Vec3[], Symbol[], Vec3[], Bool[]
    for l in content_lines(io)
        w = split(l)
        key = lowercase(w[1])
        if key == "lattice_vector"
            push!(vecs, Vec3(parse.(Float64, w[2:4])))
        elseif key in ("atom", "atom_frac")
            push!(cart, Vec3(parse.(Float64, w[2:4])))
            push!(species, species_token(w[5]))
            push!(isfrac, key == "atom_frac")
        end
    end
    isempty(species) && throw(ArgumentError("no 'atom' lines"))
    if isempty(vecs)
        any(isfrac) && throw(ArgumentError("atom_frac needs lattice_vector lines"))
        return boxed_molecule(species, cart; title)
    end
    length(vecs) == 3 || throw(ArgumentError("expected 3 lattice_vector lines, got $(length(vecs))"))
    L = Lattice(vecs...)
    return Crystal(L, species, [f ? r : cart_to_frac(L, r) for (r, f) in zip(cart, isfrac)]; title)
end

"""
    write_aims(path or io, crystal)

FHI-aims `geometry.in`: lattice vectors and fractional positions (`atom_frac`);
Cartesian `atom` lines without lattice vectors for a non-periodic structure.
"""
write_aims(path::AbstractString, c::Crystal) = open(io -> write_aims(io, c), path, "w")

function write_aims(io::IO, c::Crystal)
    println(io, "# ", isempty(c.title) ? formula(c) : c.title, " — written by AtomBond")
    if !any(c.pbc)
        foreach(((s, r),) -> @printf(io, "atom %16.10f %16.10f %16.10f %s\n", r..., s),
                zip(c.species, cartesian_positions(c)))
        return
    end
    foreach(v -> @printf(io, "lattice_vector %16.10f %16.10f %16.10f\n", v...), cell_vectors(c.lattice))
    foreach(((s, f),) -> @printf(io, "atom_frac %14.10f %14.10f %14.10f %s\n", f..., s), zip(c.species, c.frac))
end

# --- SIESTA --------------------------------------------------------------------

# fdf labels ignore case, '-', '_' and '.'.
fdf_key(s) = lowercase(replace(s, r"[-_.]" => ""))

const FDF_UNITS = Dict("ang" => 1.0, "angstrom" => 1.0, "bohr" => BOHR, "au" => BOHR, "nm" => 10.0)

"""
    read_siesta(path or io) -> Crystal

SIESTA `.fdf`: ChemicalSpeciesLabel, LatticeConstant with LatticeVectors or
LatticeParameters, AtomicCoordinatesFormat (Ang, Bohr, ScaledCartesian,
Fractional / ScaledByLatticeVectors, NotScaledCartesianAng/Bohr) and
AtomicCoordinatesAndAtomicSpecies. Without a lattice the atoms are boxed.
"""
read_siesta(path::AbstractString) = open(io -> read_siesta(io; title=splitext(basename(path))[1]), path)

function read_siesta(io::IO; title::AbstractString="")
    labels, blocks = Dict{String,Vector{String}}(), Dict{String,Vector{String}}()
    current = nothing
    for l in content_lines(io)
        w = split(l)
        if lowercase(w[1]) == "%block"
            current = fdf_key(w[2]); blocks[current] = String[]
        elseif lowercase(w[1]) == "%endblock"
            current = nothing
        elseif current !== nothing
            push!(blocks[current], l)
        else
            labels[fdf_key(w[1])] = w[2:end]
        end
    end
    haskey(blocks, "chemicalspecieslabel") || throw(ArgumentError("no ChemicalSpeciesLabel block"))
    spc = Dict(parse(Int, split(l)[1]) => species_token(split(l)[2]) for l in blocks["chemicalspecieslabel"])
    lc = get(labels, "latticeconstant", ["1.0", "Ang"])
    alat = parse(Float64, lc[1]) * get(FDF_UNITS, lowercase(get(lc, 2, "ang")), 1.0)
    L = if haskey(blocks, "latticevectors")
        Lattice((alat * Vec3(parse.(Float64, split(blocks["latticevectors"][i])[1:3])) for i in 1:3)...)
    elseif haskey(blocks, "latticeparameters")
        p = parse.(Float64, split(blocks["latticeparameters"][1])[1:6])
        Lattice((alat .* p[1:3])..., p[4:6]...)
    else
        nothing
    end
    haskey(blocks, "atomiccoordinatesandatomicspecies") || throw(ArgumentError("no AtomicCoordinatesAndAtomicSpecies block"))
    rows = [split(l) for l in blocks["atomiccoordinatesandatomicspecies"]]
    species = [spc[parse(Int, w[4])] for w in rows]
    xyz = [Vec3(parse.(Float64, w[1:3])) for w in rows]
    fmt = lowercase(get(labels, "atomiccoordinatesformat", ["Bohr"])[1])
    if fmt in ("fractional", "scaledbylatticevectors")
        L === nothing && throw(ArgumentError("fractional coordinates need a lattice"))
        return Crystal(L, species, xyz; title)
    end
    s = fmt in ("ang", "notscaledcartesianang") ? 1.0 : fmt == "scaledcartesian" ? alat : BOHR
    r = [s * x for x in xyz]
    return L === nothing ? boxed_molecule(species, r; title) : crystal_from_cartesian(L, species, r; title)
end

# --- ABINIT --------------------------------------------------------------------

"""Values after ABINIT variable `name` (no dataset suffix), with "n*x" expanded; and a unit word if given."""
function abinit_values(tokens, name)
    k = findfirst(==(name), tokens)
    k === nothing && return nothing, nothing
    vals, unit = String[], nothing
    for t in tokens[k+1:end]
        if (m = match(r"^(\d+)\*(.+)$", t)) !== nothing
            append!(vals, fill(m[2], parse(Int, m[1])))
        elseif tryparse(Float64, replace(lowercase(t), 'd' => 'e')) !== nothing || occursin(r"^[\d.+-]+/[\d.]+$", t) ||
               lowercase(t) in ("sqrt(3)", "sqrt(2)")
            push!(vals, t)
        elseif lowercase(t) in ("angstrom", "angstroms", "angstr", "bohr", "bohrs", "au", "nm")
            unit = lowercase(t)
            break
        else
            break
        end
    end
    return vals, unit
end

function abinit_number(t)
    t = lowercase(t)
    t == "sqrt(3)" && return sqrt(3.0)
    t == "sqrt(2)" && return sqrt(2.0)
    occursin('/', t) && return (/)(parse.(Float64, split(t, '/'))...)
    return parse(Float64, replace(t, 'd' => 'e'))
end

abinit_scale(unit) = unit === nothing || startswith(unit, "bohr") || unit == "au" ? BOHR : unit == "nm" ? 10.0 : 1.0

"""
    read_abinit(path or io) -> Crystal

ABINIT input: acell (Bohr unless a unit is given) with rprim or angdeg,
natom, typat, znucl, and xred, xcart or xangst. Dataset-specific variables
(acell2, ...) are ignored.
"""
read_abinit(path::AbstractString) = open(io -> read_abinit(io; title=splitext(basename(path))[1]), path)

function read_abinit(io::IO; title::AbstractString="")
    # Comments start with # or !.
    tokens = lowercase.(split(join((split(split(l, '#')[1], '!')[1] for l in eachline(io)), " ")))
    get1(name, default=nothing) = (v = abinit_values(tokens, name)[1]; v === nothing ? default : v)
    acell_v, acell_u = abinit_values(tokens, "acell")
    acell = acell_v === nothing ? [1.0, 1.0, 1.0] : abinit_scale(acell_u) .* abinit_number.(acell_v[1:3])
    L = if (ad = get1("angdeg")) !== nothing
        Lattice(acell..., abinit_number.(ad[1:3])...)
    else
        rp = get1("rprim", ["1", "0", "0", "0", "1", "0", "0", "0", "1"])
        Lattice((acell[i] * Vec3(abinit_number.(rp[3i-2:3i])) for i in 1:3)...)
    end
    natom = parse(Int, get1("natom", ["1"])[1])
    typat = parse.(Int, get1("typat", fill("1", natom))[1:natom])
    znucl = get1("znucl")
    znucl === nothing && throw(ArgumentError("no znucl"))
    species = [element(round(Int, abinit_number(znucl[t]))).symbol for t in typat]
    if (x = get1("xred")) !== nothing
        return Crystal(L, species, [Vec3(abinit_number.(x[3i-2:3i])) for i in 1:natom]; title)
    end
    x, u = abinit_values(tokens, "xcart")
    s = abinit_scale(u)
    if x === nothing
        x, s = get1("xangst"), 1.0
        x === nothing && throw(ArgumentError("no xred, xcart or xangst"))
    end
    return crystal_from_cartesian(L, species, [s * Vec3(abinit_number.(x[3i-2:3i])) for i in 1:natom]; title)
end
