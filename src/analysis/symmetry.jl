# Space-group analysis and standard cells via spglib.

"""
Crystal → spglib cell. Atom kinds are numbered by first appearance (not Z),
so unknown labels stay distinct. A kind is (species, magnetic moment,
charge), each rounded to 0.01: atoms with opposite moments (antiferromagnets)
or different charges are not symmetry-equivalent, and primitive /
conventional cells keep them apart.
"""
function spglib_cell(c::Crystal)
    key(i) = (c.species[i], round(c.magmoms[i]; digits=2), round(c.charges[i]; digits=2))
    kinds = unique(key(i) for i in 1:natoms(c))
    types = [findfirst(isequal(key(i)), kinds) for i in 1:natoms(c)]
    return Spglib.SpglibCell(Matrix(c.lattice.matrix), [Vector(f) for f in c.frac], types), kinds
end

function from_spglib(cell, kinds, title::AbstractString)
    k = kinds[cell.atoms]
    Crystal(Lattice(Matrix(cell.lattice)), first.(k), [Vec3(p) for p in cell.positions]; title,
            magmoms=[x[2] for x in k], charges=[x[3] for x in k])
end

"""
    symmetry_dataset(c; symprec=1e-3) -> Spglib.Dataset or nothing

Full spglib dataset; `nothing` for structures that are not periodic in all
three directions.
"""
function symmetry_dataset(c::Crystal; symprec::Real=1e-3)
    all(c.pbc) || return nothing
    cell, _ = spglib_cell(c)
    return Spglib.get_dataset(cell, symprec)
end

"""
    spacegroup(c; symprec=1e-3) -> (number, symbol) or nothing

e.g. `(227, "Fd-3m")` for silicon. `symprec` is the position tolerance in Å.
"""
function spacegroup(c::Crystal; symprec::Real=1e-3)
    d = symmetry_dataset(c; symprec)
    d === nothing && return nothing
    return (Int(d.spacegroup_number), String(d.international_symbol))
end

"""
    primitive_cell(c; symprec=1e-3) -> Crystal

Smallest cell with the same lattice (spglib's standardized primitive cell).
"""
function primitive_cell(c::Crystal; symprec::Real=1e-3)
    all(c.pbc) || throw(ArgumentError("primitive_cell needs a structure periodic in all directions"))
    cell, spc = spglib_cell(c)
    return from_spglib(Spglib.standardize_cell(cell, symprec; to_primitive=true), spc, c.title)
end

"""
    conventional_cell(c; symprec=1e-3) -> Crystal

Standard conventional cell for the space group (e.g. the 8-atom cubic cell
of silicon), with symmetry-idealized positions.
"""
function conventional_cell(c::Crystal; symprec::Real=1e-3)
    all(c.pbc) || throw(ArgumentError("conventional_cell needs a structure periodic in all directions"))
    cell, spc = spglib_cell(c)
    return from_spglib(Spglib.standardize_cell(cell, symprec), spc, c.title)
end
