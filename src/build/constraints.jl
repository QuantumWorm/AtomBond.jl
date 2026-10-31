# Constraints for geometry optimization, and a periodic box for molecules.
#
# Atom constraints fix fractional coordinates of atoms (VASP selective
# dynamics). Cell constraints follow CASTEP's CELL_CONSTRAINTS: six numbers for
# a b c α β γ; 0 = fixed, equal positive numbers = varied together. Both are
# stored in the Crystal and written by the VASP, CASTEP and QE exporters, as
# exactly as each code allows (see `cell_constraints_vasp`, `cell_constraints_qe`).

"""
    fix_atoms(c, indices, axes=(true, true, true)) -> Crystal

Fix the fractional coordinates `axes` (x, y, z = along a, b, c) of the given
atoms; the other coordinates of those atoms become free. `axes =
(false, false, false)` frees them.
"""
function fix_atoms(c::Crystal, idx, axes=(true, true, true))
    idx = check_indices(c, idx)
    fixed = copy(c.fixed)
    fixed[idx] .= Ref(Tuple(Bool.(axes)))
    return with_constraints(c; fixed)
end

"""Free all atom coordinates."""
free_atoms(c::Crystal) = with_constraints(c; fixed=fill(FREE_ATOM, natoms(c)))

"""
Cell constraint presets: (name, (a, b, c, α, β, γ)); 0 = fixed, equal
numbers = varied together.
"""
const CELL_CONSTRAINT_PRESETS = [
    ("All free (full cell relaxation)", (1, 2, 3, 4, 5, 6)),
    ("Fixed cell (relax atoms only)", (0, 0, 0, 0, 0, 0)),
    ("Angles fixed (relax a, b, c)", (1, 2, 3, 0, 0, 0)),
    ("Only c relaxes (a, b fixed: strained layer)", (0, 0, 1, 0, 0, 0)),
    ("a = b and c relax, angles fixed (tetragonal, hexagonal)", (1, 1, 2, 0, 0, 0)),
    ("Volume only: a b c scale together, angles fixed", (1, 1, 1, 0, 0, 0)),
]

"""
Renumber cell constraints so tied groups get 1, 2, … in order of first
appearance: (5, 5, 9, 0, 0, 0) → (1, 1, 2, 0, 0, 0). CASTEP requires the
angle numbers to differ from the length numbers, which this keeps.
"""
function normalize_cell_constraints(cc)
    ids = Dict{Int,Int}()
    return ntuple(6) do k
        cc[k] == 0 ? 0 : get!(ids, cc[k], length(ids) + 1)
    end
end

"""Lengths and angles tied: a length and an angle cannot share a number."""
function check_cell_constraints(cc)
    length(cc) == 6 && all(>=(0), cc) || throw(ArgumentError("give six numbers ≥ 0 for a b c α β γ"))
    isempty(intersect(filter(>(0), collect(cc[1:3])), filter(>(0), collect(cc[4:6])))) ||
        throw(ArgumentError("a length and an angle cannot be tied (use different numbers)"))
    return normalize_cell_constraints(cc)
end

"""Words for cell constraints: "a, b fixed; c free; angles fixed"."""
function describe_cell_constraints(cc)
    cc == (0, 0, 0, 0, 0, 0) && return "cell fixed"
    names = ("a", "b", "c", "α", "β", "γ")
    parts = String[]
    for (group, label) in (((1, 2, 3), "lengths"), ((4, 5, 6), "angles"))
        vals = [cc[k] for k in group]
        if all(==(0), vals)
            push!(parts, "$label fixed")
        elseif length(unique(vals)) == 3
            push!(parts, "$label free")
        else
            fixed = [names[k] for k in group if cc[k] == 0]
            isempty(fixed) || push!(parts, join(fixed, ", ") * " fixed")
            for n in unique(filter(>(0), vals))
                ks = [names[k] for k in group if cc[k] == n]
                push!(parts, length(ks) == 1 ? "$(ks[1]) free" : join(ks, " = ") * " together")
            end
        end
    end
    return join(parts, "; ")
end

isorthogonal(L::Lattice) = all(isapprox.(lattice_parameters(L)[4:6], 90; atol=1e-6))

"""
    cell_constraints_vasp(c) -> (ISIF or nothing, note)

VASP's ISIF for the cell constraints: 2 = cell fixed, 3 = all free, 7 =
volume only. Other constraints need a patched VASP; `nothing` then.
"""
function cell_constraints_vasp(c::Crystal)
    cc = c.cell_constraints
    cc == (0, 0, 0, 0, 0, 0) && return 2, "ISIF = 2: cell fixed"
    cc == FREE_CELL && return 3, "ISIF = 3: cell and atoms relax"
    cc == (1, 1, 1, 0, 0, 0) && return 7, "ISIF = 7: volume only"
    return nothing, "VASP cannot fix single lattice parameters without a patch (ISIF = 3 relaxes all)"
end

"""
    cell_constraints_qe(c) -> (calculation, cell_dofree or nothing, note)

pw.x settings for the constraints: 'relax' for a fixed cell, else
'vc-relax' with the matching `cell_dofree` when one exists.
"""
function cell_constraints_qe(c::Crystal)
    cc = c.cell_constraints
    cc == (0, 0, 0, 0, 0, 0) && return "relax", nothing, "calculation = 'relax' (cell fixed)"
    dofree = cc == FREE_CELL ? "all" :
             cc == (1, 1, 1, 0, 0, 0) ? "volume" :
             cc == (1, 2, 3, 0, 0, 0) && isorthogonal(c.lattice) && is_standard_axes(c.lattice) ? "xyz" :
             cc == (0, 0, 1, 0, 0, 0) && is_standard_axes(c.lattice) && c_along_z(c.lattice) ? "z" :
             nothing
    dofree === nothing &&
        return "vc-relax", "all", "pw.x has no cell_dofree for these constraints; written as 'all' (edit by hand)"
    return "vc-relax", dofree, "calculation = 'vc-relax', cell_dofree = '$dofree'"
end

is_standard_axes(L::Lattice) = (M = L.matrix; abs(M[2, 1]) < 1e-8 && abs(M[3, 1]) < 1e-8 && abs(M[3, 2]) < 1e-8)
c_along_z(L::Lattice) = (M = L.matrix; abs(M[1, 3]) < 1e-8 && abs(M[2, 3]) < 1e-8)

"""
    periodic_box(c; vacuum=10.0, cubic=false) -> Crystal

Put a molecule (a structure that is not periodic) into a periodic box for
plane-wave codes (CASTEP, VASP, QE): an orthorhombic cell `vacuum` Å larger
than the molecule along x, y and z (`cubic`: all three equal), with the
molecule in the middle. The result is periodic, so its cell is shown and
exported; k-points along the vacuum directions are 1 (Γ).
"""
function periodic_box(c::Crystal; vacuum::Real=10.0, cubic::Bool=false)
    any(c.pbc) && throw(ArgumentError("only a molecule (non-periodic structure) can be put in a box"))
    vacuum > 0 || throw(ArgumentError("vacuum must be positive"))
    r = cartesian_positions(c)
    lo = isempty(r) ? Vec3(0, 0, 0) : reduce((a, b) -> min.(a, b), r)
    hi = isempty(r) ? Vec3(0, 0, 0) : reduce((a, b) -> max.(a, b), r)
    size = hi - lo .+ vacuum
    cubic && (size = Vec3(fill(maximum(size), 3)...))
    L = Lattice(Vec3(size[1], 0, 0), Vec3(0, size[2], 0), Vec3(0, 0, size[3]))
    shift = size / 2 - (lo + hi) / 2
    return Crystal(L, c.species, [cart_to_frac(L, p + shift) for p in r]; pbc=(true, true, true), title=c.title,
                   fixed=c.fixed, cell_constraints=c.cell_constraints, charges=c.charges,
                   magmoms=c.magmoms)
end

"""A periodic structure back to a molecule (not periodic; the cell is only a container)."""
as_molecule(c::Crystal) = Crystal(c.lattice, c.species, c.frac; pbc=(false, false, false), title=c.title,
                                  fixed=c.fixed, cell_constraints=c.cell_constraints, charges=c.charges,
                   magmoms=c.magmoms)
