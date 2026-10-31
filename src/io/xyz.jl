# XYZ and extended XYZ (Lattice="..." in the comment line).
#
# A plain XYZ file has no cell, so the molecule is placed in an orthorhombic
# box with `padding` Å of empty space on each side and pbc = (false, false, false).

"""
    read_xyz(path; padding=5.0) -> Crystal

Reads the first frame only.
"""
read_xyz(path::AbstractString; kwargs...) = open(io -> read_xyz(io; name=basename(path), kwargs...), path)

function read_xyz(io::IO; name::AbstractString="xyz", padding::Real=5.0)
    n = parse(Int, strip(readline(io)))
    comment = readline(io)
    species = Symbol[]
    positions = Vec3[]
    qcol = charge_column(comment)
    charges = qcol === nothing ? nothing : Float64[]
    mcol = charge_column(comment; names=("magmom", "magmoms", "initial_magmoms", "magnetic_moment", "magnetic_moments"))
    magmoms = mcol === nothing ? nothing : Float64[]
    for i in 1:n
        eof(io) && error("$name: expected $n atoms, found $(i-1)")
        t = split(readline(io))
        push!(species, normalize_symbol(t[1]))
        push!(positions, Vec3(parse.(Float64, t[2:4])))
        qcol === nothing || push!(charges, parse(Float64, t[qcol]))
        mcol === nothing || push!(magmoms, parse(Float64, t[mcol]))
    end

    m = match(r"Lattice\s*=\s*\"([^\"]+)\""i, comment)
    if m !== nothing
        v = parse.(Float64, split(m.captures[1]))
        length(v) == 9 || error("$name: Lattice= must have 9 numbers")
        L = Lattice(reshape(v, 3, 3))   # extxyz lists a, b, c in order → columns
        pm = match(r"pbc\s*=\s*\"([^\"]+)\""i, comment)
        pbc = pm === nothing ? (true, true, true) :
              Tuple(lowercase(t) in ("t", "true", "1") for t in split(pm.captures[1]))
        # Whatever is left after removing key=value pairs is the title.
        title = strip(replace(comment, r"\w+\s*=\s*\"[^\"]*\"" => "", r"\w+=\S+" => ""))
        c = crystal_from_cartesian(L, species, positions; pbc=pbc, title=title)
        return with_magmoms(with_charges(c, charges), magmoms)
    end

    return with_magmoms(with_charges(boxed_molecule(species, positions; padding, title=strip(comment)), charges), magmoms)
end

"""
Column of a per-atom property in an extended XYZ file (`Properties=
name:type:count:...`), by default the charge (charge, charges or
initial_charges), or `nothing`.
"""
function charge_column(comment::AbstractString; names=("charge", "charges", "initial_charges"))
    m = match(r"Properties\s*=\s*(\S+)"i, comment)
    m === nothing && return nothing
    parts = split(m[1], ':')
    col = 1
    for k in 1:3:length(parts)-2
        lowercase(parts[k]) in names && return col
        col += parse(Int, parts[k+2])
    end
    return nothing
end

"""
    boxed_molecule(species, positions; padding=5.0, title="") -> Crystal

A non-periodic structure (molecule, cluster) in an orthorhombic box with
`padding` Å of space around the atoms, for formats without a cell.
"""
function boxed_molecule(species, positions; padding::Real=5.0, title::AbstractString="")
    isempty(positions) && throw(ArgumentError("no atoms"))
    lo = reduce((a, b) -> min.(a, b), positions)
    hi = reduce((a, b) -> max.(a, b), positions)
    L = Lattice(Diagonal((hi .- lo) .+ 2padding))
    shifted = [Vec3(r .- lo .+ padding) for r in positions]
    return crystal_from_cartesian(L, species, shifted; pbc=(false, false, false), title)
end

"""
    write_xyz(path, crystal)

Write extended XYZ with the cell in the comment line.
"""
write_xyz(path::AbstractString, c::Crystal) = open(io -> write_xyz(io, c), path, "w")

function write_xyz(io::IO, c::Crystal)
    println(io, natoms(c))
    lat = join((@sprintf("%.10f", x) for x in vec(c.lattice.matrix)), " ")
    pbc = join((p ? "T" : "F" for p in c.pbc), " ")
    charged, magnetic = has_charges(c), has_magmoms(c)      # unknown values are written as 0
    println(io, "Lattice=\"$lat\" Properties=species:S:1:pos:R:3", charged ? ":charge:R:1" : "",
            magnetic ? ":magmom:R:1" : "", " pbc=\"$pbc\"")
    for (i, (s, r)) in enumerate(zip(c.species, cartesian_positions(c)))
        @printf(io, "%-3s %16.10f %16.10f %16.10f", s, r...)
        charged && @printf(io, " %10.6f", isnan(c.charges[i]) ? 0.0 : c.charges[i])
        magnetic && @printf(io, " %10.6f", isnan(c.magmoms[i]) ? 0.0 : c.magmoms[i])
        println(io)
    end
end
