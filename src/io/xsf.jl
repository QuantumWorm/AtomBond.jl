# XCrySDen XSF (also read and written by VESTA, Quantum ESPRESSO, ...) and
# Gaussian cube files (the atoms and cell; the volumetric data is skipped).

"""Element of a species token: symbol ("Si", "SI", "Si1"), labelled ("Si:1") or atomic number ("14")."""
species_token(t::AbstractString) = castep_species(t)

# Comment-free, stripped, non-empty lines.
function content_lines(io::IO; comment::AbstractString="#")
    lines = String[]
    for raw in eachline(io)
        line = strip(split(raw, comment)[1])
        isempty(line) || push!(lines, line)
    end
    return lines
end

"""
    read_xsf(path or io) -> Crystal

XSF file: CRYSTAL / SLAB / POLYMER (periodic in 3 / 2 / 1 directions, cell
from PRIMVEC, atoms from PRIMCOORD) or MOLECULE (ATOMS, boxed). Atoms may be
given by symbol or atomic number; forces after the coordinates are ignored.
For animated files (ANIMSTEPS) the last step is read.
"""
read_xsf(path::AbstractString) = open(io -> read_xsf(io; title=splitext(basename(path))[1]), path)

function read_xsf(io::IO; title::AbstractString="")
    lines = content_lines(io)
    kind = :crystal
    L = nothing
    species, positions = Symbol[], Vec3[]
    keywords = ("CRYSTAL", "SLAB", "POLYMER", "MOLECULE", "PRIMVEC", "CONVVEC", "PRIMCOORD", "CONVCOORD",
                "ATOMS", "ANIMSTEPS", "BEGIN_", "END_")
    iskey(l) = any(k -> startswith(uppercase(l), k), keywords)
    atom(w) = (push!(species, species_token(w[1])); push!(positions, Vec3(parse.(Float64, w[2:4]))))
    k = 1
    while k <= length(lines)
        key = uppercase(split(lines[k])[1])
        if key in ("CRYSTAL", "SLAB", "POLYMER", "MOLECULE")
            kind = Symbol(lowercase(key))
        elseif key == "PRIMVEC"
            L = Lattice((Vec3(parse.(Float64, split(lines[k+i])[1:3])) for i in 1:3)...)
            k += 3
        elseif key == "PRIMCOORD"
            n = parse(Int, split(lines[k+1])[1])
            empty!(species); empty!(positions)
            foreach(i -> atom(split(lines[k+1+i])), 1:n)
            k += n + 1
        elseif key == "ATOMS"
            empty!(species); empty!(positions)
            while k + 1 <= length(lines) && !iskey(lines[k+1])
                atom(split(lines[k+1]))
                k += 1
            end
        elseif key == "CONVVEC"
            k += 3
        elseif key == "CONVCOORD"
            k += parse(Int, split(lines[k+1])[1]) + 1
        elseif startswith(key, "BEGIN_BLOCK") || startswith(key, "BEGIN_DATAGRID")
            break                            # volumetric data follows
        end
        k += 1
    end
    isempty(positions) && throw(ArgumentError("no atoms found (PRIMCOORD or ATOMS)"))
    if kind === :molecule || L === nothing
        return boxed_molecule(species, positions; title)
    end
    pbc = kind === :slab ? (true, true, false) : kind === :polymer ? (true, false, false) : (true, true, true)
    return crystal_from_cartesian(L, species, positions; pbc, title)
end

"""
    write_xsf(path or io, crystal)

XSF with PRIMVEC and PRIMCOORD (Cartesian Å): CRYSTAL, SLAB or POLYMER after
the periodic directions, MOLECULE (ATOMS) for a non-periodic structure.
"""
write_xsf(path::AbstractString, c::Crystal) = open(io -> write_xsf(io, c), path, "w")

function write_xsf(io::IO, c::Crystal)
    println(io, "# ", isempty(c.title) ? formula(c) : c.title, " — written by AtomBond")
    r = cartesian_positions(c)
    if !any(c.pbc)
        println(io, "MOLECULE\nATOMS")
        foreach(i -> @printf(io, "%-3s %16.10f %16.10f %16.10f\n", c.species[i], r[i]...), eachindex(r))
        return
    end
    println(io, c.pbc == (true, true, false) ? "SLAB" : c.pbc == (true, false, false) ? "POLYMER" : "CRYSTAL")
    println(io, "PRIMVEC")
    foreach(v -> @printf(io, "  %16.10f %16.10f %16.10f\n", v...), cell_vectors(c.lattice))
    println(io, "PRIMCOORD\n", natoms(c), " 1")
    foreach(i -> @printf(io, "%-3s %16.10f %16.10f %16.10f\n", c.species[i], r[i]...), eachindex(r))
end

"""
    read_cube(path or io) -> Crystal

Gaussian cube file: the atoms and the cell spanned by the grid (Bohr, or Å
when the grid counts are negative). The volumetric data is not read.
"""
read_cube(path::AbstractString) = open(io -> read_cube(io; title=splitext(basename(path))[1]), path)

function read_cube(io::IO; title::AbstractString="")
    c1 = strip(readline(io)); readline(io)
    w = split(readline(io))
    nat = parse(Int, w[1])
    origin = Vec3(parse.(Float64, w[2:4]))
    axes = map(1:3) do _
        a = split(readline(io))
        n = parse(Int, a[1])
        (n, Vec3(parse.(Float64, a[2:4])))
    end
    unit = axes[1][1] < 0 ? 1.0 : BOHR
    L = Lattice((unit * abs(n) * v for (n, v) in axes)...)
    species, positions = Symbol[], Vec3[]
    for _ in 1:abs(nat)
        a = split(readline(io))
        push!(species, species_token(a[1]))
        push!(positions, unit * (Vec3(parse.(Float64, a[3:5])) - origin))
    end
    return crystal_from_cartesian(L, species, positions; title=isempty(c1) ? title : c1)
end
