# PDB and LAMMPS data files (read and write).

# --- PDB -------------------------------------------------------------------------

"""
    read_pdb(path or io) -> Crystal

PDB: the cell from CRYST1 (a structure without one, or with the 1 Å
placeholder cell, is boxed as a molecule) and ATOM / HETATM records of the
first MODEL. The element comes from columns 77–78, otherwise from the atom
name; a formal charge in columns 79–80 ("2+", "1-") is kept.
"""
read_pdb(path::AbstractString) = open(io -> read_pdb(io; title=splitext(basename(path))[1]), path)

function read_pdb(io::IO; title::AbstractString="")
    cell = nothing
    species, positions, charges = Symbol[], Vec3[], Float64[]
    col(l, r) = length(l) >= last(r) ? strip(l[r]) : length(l) >= first(r) ? strip(l[first(r):end]) : ""
    for l in eachline(io)
        rec = uppercase(col(l, 1:6))
        if rec == "CRYST1"
            p = parse.(Float64, split(l[7:min(end, 54)]))
            cell = all(≈(1.0), p[1:3]) ? nothing : p          # 1 1 1: no cell (PDB convention)
        elseif rec in ("ATOM", "HETATM")
            el = col(l, 77:78)
            el = isempty(el) || !isletter(el[1]) ? col(l, 13:16) : el
            push!(species, normalize_symbol(el))
            xyz = tryparse.(Float64, (col(l, 31:38), col(l, 39:46), col(l, 47:54)))
            any(isnothing, xyz) && (xyz = parse.(Float64, split(l)[end-4:end-2]))   # loosely formatted
            push!(positions, Vec3(xyz...))
            q = match(r"^(\d)([+-])$", col(l, 79:80))
            push!(charges, q === nothing ? NaN : (q[2] == "-" ? -1 : 1) * parse(Int, q[1]))
        elseif rec in ("ENDMDL", "END")
            isempty(positions) || break
        elseif startswith(rec, "TITLE") || startswith(rec, "COMPND")
            isempty(title) && (title = col(l, 11:80))
        end
    end
    isempty(positions) && throw(ArgumentError("no ATOM or HETATM records"))
    c = cell === nothing ? boxed_molecule(species, positions; title) :
        crystal_from_cartesian(Lattice(cell...), species, positions; title)
    return any(!isnan, charges) ? with_charges(c, charges) : c
end

"""
    write_pdb(path or io, crystal)

PDB with CRYST1 (space group P 1) and HETATM records. Positions are Cartesian
in the standard orientation (a along x, b in the xy plane), as PDB requires.
"""
write_pdb(path::AbstractString, c::Crystal) = open(io -> write_pdb(io, c), path, "w")

function write_pdb(io::IO, c::Crystal)
    println(io, rpad("TITLE     " * (isempty(c.title) ? formula(c) : c.title), 80)[1:80])
    a, b, cc, α, β, γ = lattice_parameters(c.lattice)
    any(c.pbc) && @printf(io, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n", a, b, cc, α, β, γ, "P 1", 1)
    L = Lattice(a, b, cc, α, β, γ)
    count = Dict{Symbol,Int}()
    for (i, (s, f)) in enumerate(zip(c.species, c.frac))
        r = frac_to_cart(L, f)
        n = count[s] = get(count, s, 0) + 1
        name = first(string(s, n), 4)
        @printf(io, "HETATM%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
                i % 100000, name, "MOL", "A", 1, r..., 1.0, 0.0, uppercase(string(s)))
    end
    println(io, "END")
end

# --- LAMMPS ------------------------------------------------------------------------

"""Element with atomic mass closest to `m` (within 0.5 u)."""
function element_by_mass(m::Real)
    best = argmin(e -> abs(e.mass - m), collect(values(ELEMENTS)))
    abs(best.mass - m) < 0.5 || throw(ArgumentError("no element with mass $m"))
    return best.symbol
end

"""
    read_lammps_data(path or io; style=nothing) -> Crystal

LAMMPS data file: box (with xy xz yz tilts), Masses and Atoms (atomic,
charge or full style; taken from the "Atoms # style" comment or the number
of columns). Each type's element is the element named in the Masses comment
(`1 28.0855 # Si`), otherwise the one with the closest mass.
"""
read_lammps_data(path::AbstractString; kw...) =
    open(io -> read_lammps_data(io; title=splitext(basename(path))[1], kw...), path)

function read_lammps_data(io::IO; title::AbstractString="", style=nothing)
    lines = readlines(io)
    isempty(title) && (title = strip(lines[1]))
    natoms = 0
    lo, hi, tilt = zeros(3), zeros(3), zeros(3)
    types = Dict{Int,Symbol}()
    rows = Vector{SubString{String}}[]
    k = 2
    section(name) = occursin(Regex("^\\s*$name\\b"), lines[k])
    while k <= length(lines)
        l = split(lines[k], '#')[1]
        w = split(l)
        if occursin(r"\batoms\s*$", l)
            natoms = parse(Int, w[1])
        elseif occursin(r"xlo\s+xhi", l)
            lo[1], hi[1] = parse.(Float64, w[1:2])
        elseif occursin(r"ylo\s+yhi", l)
            lo[2], hi[2] = parse.(Float64, w[1:2])
        elseif occursin(r"zlo\s+zhi", l)
            lo[3], hi[3] = parse.(Float64, w[1:2])
        elseif occursin(r"xy\s+xz\s+yz", l)
            tilt .= parse.(Float64, w[1:3])
        elseif section("Masses")
            k += 1
            while k <= length(lines) && (isempty(strip(lines[k])) || tryparse(Int, split(lines[k])[1]) !== nothing)
                if !isempty(strip(lines[k]))
                    parts = split(lines[k], '#')
                    t, m = split(parts[1])[1:2]
                    word = length(parts) > 1 && !isempty(split(parts[2])) ? first(split(parts[2])) : ""
                    el = normalize_symbol(word)
                    named = haskey(ELEMENTS, el) && lowercase(string(el)) == lowercase(word)   # "Si", not "Si_bulk"
                    types[parse(Int, t)] = named ? el : element_by_mass(parse(Float64, m))
                end
                k += 1
            end
            continue
        elseif section("Atoms")
            m = match(r"#\s*(\w+)", lines[k])
            style === nothing && m !== nothing && (style = m[1])
            k += 1
            while k <= length(lines) && length(rows) < natoms
                isempty(strip(lines[k])) || push!(rows, split(split(lines[k], '#')[1]))
                k += 1
            end
            continue
        end
        k += 1
    end
    length(rows) == natoms > 0 || throw(ArgumentError("expected $natoms atoms, found $(length(rows))"))
    n = length(rows[1])
    style = something(style, n in (5, 8) ? "atomic" : n in (6, 9) ? "charge" : n in (7, 10) ? "full" : "atomic")
    tcol, xcol = style == "full" ? (3, 5) : style == "charge" ? (2, 4) : style == "atomic" ? (2, 3) :
                 throw(ArgumentError("atom style '$style' is not supported (atomic, charge, full)"))
    sort!(rows; by=r -> parse(Int, r[1]))
    L = Lattice(Vec3(hi[1] - lo[1], 0, 0), Vec3(tilt[1], hi[2] - lo[2], 0), Vec3(tilt[2], tilt[3], hi[3] - lo[3]))
    species = map(rows) do r
        t = parse(Int, r[tcol])
        haskey(types, t) || throw(ArgumentError("atom type $t has no Masses entry"))
        types[t]
    end
    r = [Vec3(parse.(Float64, w[xcol:xcol+2])) - Vec3(lo) for w in rows]
    # charge and full styles have the charge q just before x.
    charges = style == "atomic" ? nothing : [parse(Float64, w[xcol-1]) for w in rows]
    return with_charges(crystal_from_cartesian(L, species, r; title), charges)
end

"""
    write_lammps_data(path or io, crystal)

LAMMPS data file, atomic style (charge style when the structure has
charges; unknown ones are written as 0): the cell in LAMMPS' restricted
triclinic form (a along x, b in the xy plane; tilts written when non-zero),
one atom type per element with the element named in the Masses comment.
"""
write_lammps_data(path::AbstractString, c::Crystal) = open(io -> write_lammps_data(io, c), path, "w")

function write_lammps_data(io::IO, c::Crystal)
    spc = species_list(c)
    L = Lattice(lattice_parameters(c.lattice)...)
    a, b, cc = cell_vectors(L)
    println(io, isempty(c.title) ? formula(c) : c.title, " — written by AtomBond\n")
    println(io, natoms(c), " atoms\n", length(spc), " atom types\n")
    @printf(io, "0.0 %.10f xlo xhi\n0.0 %.10f ylo yhi\n0.0 %.10f zlo zhi\n", a[1], b[2], cc[3])
    any(abs(x) > 1e-10 for x in (b[1], cc[1], cc[2])) && @printf(io, "%.10f %.10f %.10f xy xz yz\n", b[1], cc[1], cc[2])
    println(io, "\nMasses\n")
    foreach(((t, s),) -> @printf(io, "%d %.4f # %s\n", t, atomic_mass(s), s), enumerate(spc))
    charged = has_charges(c)
    println(io, "\nAtoms # ", charged ? "charge" : "atomic", "\n")
    for (i, (s, f)) in enumerate(zip(c.species, c.frac))
        @printf(io, "%d %d ", i, findfirst(==(s), spc))
        charged && @printf(io, "%.6f ", isnan(c.charges[i]) ? 0.0 : c.charges[i])
        @printf(io, "%.10f %.10f %.10f\n", frac_to_cart(L, f)...)
    end
end
