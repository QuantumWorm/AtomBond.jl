# Materials Studio (BIOVIA) .car archives and .xsd documents, and SHELX
# .res / .ins files (as written by AIRSS and many structure tools).

# --- Materials Studio .car ------------------------------------------------------

"""
    read_car(path or io) -> Crystal

Materials Studio / Insight `.car` archive (version 3): the PBC line gives
a b c α β γ, and atom lines give Cartesian Å positions with the element in
the 8th column. Positions are taken in Materials Studio's default
orientation (a along x, b in the xy plane). PBC=OFF files are boxed.
"""
read_car(path::AbstractString) = open(io -> read_car(io; title=splitext(basename(path))[1]), path)

function read_car(io::IO; title::AbstractString="")
    lines = readlines(io)
    periodic = any(l -> occursin(r"^\s*PBC\s*=\s*ON"i, l), lines)
    cell = nothing
    species, positions = Symbol[], Vec3[]
    k = findfirst(l -> occursin(r"^\s*PBC\s*=", l), lines)
    k === nothing && throw(ArgumentError("no PBC= line (not a .car file?)"))
    # The line after PBC= is the title (unless it is already the !DATE line).
    t = strip(get(lines, k + 1, ""))
    isempty(t) || startswith(t, "!") || startswith(uppercase(t), "PBC") ||
        startswith(t, "Materials Studio Generated") || (title = t)
    for l in lines[k+1:end]
        w = split(l)
        isempty(w) && continue
        if uppercase(w[1]) == "PBC" && length(w) >= 7
            cell = parse.(Float64, w[2:7])
        elseif lowercase(w[1]) == "end"
            isempty(positions) || break
        elseif length(w) >= 8 && all(t -> tryparse(Float64, t) !== nothing, w[2:4])
            push!(species, normalize_symbol(w[8]))
            push!(positions, Vec3(parse.(Float64, w[2:4])))
        end
    end
    isempty(positions) && throw(ArgumentError("no atoms"))
    (periodic && cell !== nothing) || return boxed_molecule(species, positions; title)
    return crystal_from_cartesian(Lattice(cell...), species, positions; title)
end

"""
    write_car(path or io, crystal)

Materials Studio `.car` archive (PBC=ON with a b c α β γ and P1, or PBC=OFF),
Cartesian positions in the standard orientation.
"""
write_car(path::AbstractString, c::Crystal) = open(io -> write_car(io, c), path, "w")

function write_car(io::IO, c::Crystal)
    periodic = any(c.pbc)
    println(io, "!BIOSYM archive 3\nPBC=", periodic ? "ON" : "OFF")
    println(io, isempty(c.title) ? formula(c) : c.title)
    println(io, "!DATE ", Dates.format(Dates.now(), "e u d HH:MM:SS yyyy"))
    p = lattice_parameters(c.lattice)
    periodic && @printf(io, "PBC%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f (P1)\n", p...)
    L = Lattice(p...)
    count = Dict{Symbol,Int}()
    for (s, f) in zip(c.species, c.frac)
        r = periodic ? frac_to_cart(L, f) : frac_to_cart(c.lattice, f)
        n = count[s] = get(count, s, 0) + 1
        @printf(io, "%-5s %14.9f %14.9f %14.9f XXXX 1      xx      %-2s  0.000\n", string(s, n), r..., s)
    end
    println(io, "end\nend")
end

# --- Materials Studio .xsd -----------------------------------------------------

xml_attr(tag, name) = (m = match(Regex("\\b$name=\"([^\"]*)\""), tag)) === nothing ? nothing : m[1]

"""
    read_xsd(path or io) -> Crystal

Materials Studio `.xsd` document: the lattice from SpaceGroup (AVector,
BVector, CVector) and the Atom3d elements (fractional XYZ in a periodic
document, Cartesian Å otherwise). Image atoms (ImageOf) are skipped; for a
mixed-occupancy site the first component is used.
"""
read_xsd(path::AbstractString) = open(io -> read_xsd(io; title=splitext(basename(path))[1]), path)

function read_xsd(io::IO; title::AbstractString="")
    xml = read(io, String)
    vec3(s) = Vec3(parse.(Float64, split(s, ','))[1:3])
    sg = match(r"<SpaceGroup\b[^>]*>", xml)
    L = sg === nothing ? nothing :
        Lattice(vec3(xml_attr(sg.match, "AVector")), vec3(xml_attr(sg.match, "BVector")), vec3(xml_attr(sg.match, "CVector")))
    species, xyz = Symbol[], Vec3[]
    for m in eachmatch(r"<Atom3d\b[^>]*>", xml)
        tag = m.match
        xml_attr(tag, "ImageOf") === nothing || continue
        p, comp = xml_attr(tag, "XYZ"), xml_attr(tag, "Components")
        (p === nothing || comp === nothing) && continue
        push!(species, normalize_symbol(first(split(comp, ','))))
        push!(xyz, vec3(p))
    end
    isempty(xyz) && throw(ArgumentError("no Atom3d elements with XYZ and Components"))
    L === nothing && return boxed_molecule(species, xyz; title)
    return Crystal(L, species, xyz; title)
end

# --- SHELX .res / .ins ---------------------------------------------------------

const SHELX_INSTRUCTIONS = Set(split("""
    TITL CELL ZERR LATT SYMM SFAC DISP UNIT FVAR REM END HKLF L.S. CGLS PLAN FMAP LIST WGHT
    BOND CONF ACTA MORE TEMP SIZE OMIT SHEL BASF EXTI SWAT AFIX HFIX PART RESI DFIX DANG SADI
    SAME SIMU DELU ISOR RIGU FLAT CHIV EADP EXYZ ANIS BLOC DAMP STIR TWIN MERG SPEC FREE HTAB
    MPLA RTAB EQIV CONN SUMP WPDB""", r"\s+"; keepempty=false))

# Centring translations of SHELX LATT |n|.
const SHELX_CENTRING = Dict(1 => [Vec3(0, 0, 0)], 2 => [Vec3(0, 0, 0), Vec3(0.5, 0.5, 0.5)],
    3 => [Vec3(0, 0, 0), Vec3(2 / 3, 1 / 3, 1 / 3), Vec3(1 / 3, 2 / 3, 2 / 3)],
    4 => [Vec3(0, 0, 0), Vec3(0, 0.5, 0.5), Vec3(0.5, 0, 0.5), Vec3(0.5, 0.5, 0)],
    5 => [Vec3(0, 0, 0), Vec3(0, 0.5, 0.5)], 6 => [Vec3(0, 0, 0), Vec3(0.5, 0, 0.5)],
    7 => [Vec3(0, 0, 0), Vec3(0.5, 0.5, 0)])

"""
    read_res(path or io) -> Crystal

SHELX `.res` / `.ins` (AIRSS output included): CELL, LATT (centring, and
inversion when positive), SYMM operations, SFAC element list and atom lines
(name, SFAC number, fractional x y z). Symmetry-equivalent atoms are generated.
"""
read_res(path::AbstractString) = open(io -> read_res(io; title=splitext(basename(path))[1]), path)

function read_res(io::IO; title::AbstractString="")
    cell, latt, ops, sfac = nothing, 1, SymOp[], Symbol[]
    sites = Tuple{Symbol,Vec3}[]
    for raw in eachline(io)
        l = strip(raw)
        isempty(l) && continue
        w = split(l)
        key = uppercase(w[1])
        if key == "TITL"
            length(w) > 1 && (title = strip(l[5:end]))
        elseif key == "CELL"
            cell = parse.(Float64, w[3:8])           # wavelength first
        elseif key == "LATT"
            latt = parse(Int, w[2])
        elseif key == "SYMM"
            push!(ops, parse_symop(l[5:end]))
        elseif key == "SFAC"
            # Element names only (the long form adds scattering-factor numbers).
            append!(sfac, normalize_symbol.(filter(t -> occursin(r"^[A-Za-z]{1,2}$", t), w[2:end])))
        elseif key in ("END", "HKLF")
            break
        elseif key ∉ SHELX_INSTRUCTIONS && length(w) >= 5 && tryparse(Int, w[2]) !== nothing &&
               all(t -> tryparse(Float64, t) !== nothing, w[3:5])
            t = parse(Int, w[2])
            1 <= t <= length(sfac) || throw(ArgumentError("atom $(w[1]): SFAC number $t out of range"))
            # SHELX adds 10 to fixed coordinates (x = 10.25 means fixed 0.25).
            f = Vec3((x -> abs(x) >= 5 ? x - 10sign(x) : x).(parse.(Float64, w[3:5])))
            push!(sites, (sfac[t], f))
        end
    end
    cell === nothing && throw(ArgumentError("no CELL line"))
    isempty(sites) && throw(ArgumentError("no atoms"))
    L = Lattice(cell...)
    I3 = SMatrix{3,3,Int}(1, 0, 0, 0, 1, 0, 0, 0, 1)
    base = [(I3, Vec3(0, 0, 0)); ops]
    latt > 0 && (base = [base; [(-W, -w) for (W, w) in base]])
    all_ops = [(W, w + t) for (W, w) in base for t in SHELX_CENTRING[abs(latt)]]
    species, frac = Symbol[], Vec3[]
    for (s, f) in sites
        for g in orbit(L, f, all_ops)
            push!(species, s); push!(frac, g)
        end
    end
    return Crystal(L, species, frac; title)
end

"""
    write_res(path or io, crystal)

SHELX `.res` in P1 (LATT -1, no SYMM), as AIRSS writes it.
"""
write_res(path::AbstractString, c::Crystal) = open(io -> write_res(io, c), path, "w")

function write_res(io::IO, c::Crystal)
    spc = species_list(c)
    println(io, "TITL ", isempty(c.title) ? formula(c) : replace(c.title, r"\s+" => "_"))
    @printf(io, "CELL 1.54180 %.6f %.6f %.6f %.4f %.4f %.4f\n", lattice_parameters(c.lattice)...)
    println(io, "LATT -1\nSFAC ", join(spc, " "))
    count = Dict{Symbol,Int}()
    for (s, f) in zip(c.species, c.frac)
        n = count[s] = get(count, s, 0) + 1
        @printf(io, "%-6s %2d %12.8f %12.8f %12.8f 1.0\n", string(s, n), findfirst(==(s), spc), f...)
    end
    println(io, "END")
end
