# Crystal: a lattice plus atoms given in fractional coordinates.
#
# `pbc` marks which directions are periodic. A molecule read from a plain XYZ
# file gets a padded box with pbc = (false, false, false); a slab would be
# (true, true, false).

#
# Constraints for geometry optimization travel with the structure:
# - `fixed[i]` = (x, y, z): which fractional coordinates of atom i are fixed
#   (VASP selective dynamics "F"); all free by default.
# - `cell_constraints` = (a, b, c, α, β, γ) as in CASTEP's CELL_CONSTRAINTS:
#   0 = fixed, equal positive numbers = tied together (a = b), different
#   numbers = independent. Default (1, 2, 3, 4, 5, 6): all free.
#
# `charges` are atomic charges in e when a file provides them (LAMMPS data,
# extended XYZ, PDB, CASTEP Mulliken analysis); NaN = not known.
# `magmoms` are collinear magnetic moments in μB (CASTEP SPIN= and Mulliken
# spins, extended XYZ); NaN = not known.

const FREE_CELL = (1, 2, 3, 4, 5, 6)
const FREE_ATOM = (false, false, false)

struct Crystal
    lattice::Lattice
    species::Vector{Symbol}
    frac::Vector{Vec3}
    pbc::NTuple{3,Bool}
    title::String
    fixed::Vector{NTuple{3,Bool}}
    cell_constraints::NTuple{6,Int}
    charges::Vector{Float64}
    magmoms::Vector{Float64}
    function Crystal(lattice::Lattice, species::AbstractVector, frac::AbstractVector;
                     pbc=(true, true, true), title::AbstractString="", fixed=nothing,
                     cell_constraints=FREE_CELL, charges=nothing, magmoms=nothing)
        length(species) == length(frac) ||
            throw(ArgumentError("got $(length(species)) species but $(length(frac)) positions"))
        syms = Symbol[s isa Symbol ? s : normalize_symbol(string(s)) for s in species]
        fx = fixed === nothing ? fill(FREE_ATOM, length(syms)) : NTuple{3,Bool}[Tuple(Bool.(f)) for f in fixed]
        length(fx) == length(syms) || throw(ArgumentError("got $(length(fx)) constraints for $(length(syms)) atoms"))
        cc = Tuple(Int.(cell_constraints))
        (length(cc) == 6 && all(>=(0), cc)) ||
            throw(ArgumentError("cell constraints are six numbers ≥ 0 (a b c α β γ), got $cell_constraints"))
        q = charges === nothing ? fill(NaN, length(syms)) : Float64[something(x, NaN) for x in charges]
        length(q) == length(syms) || throw(ArgumentError("got $(length(q)) charges for $(length(syms)) atoms"))
        m = magmoms === nothing ? fill(NaN, length(syms)) : Float64[something(x, NaN) for x in magmoms]
        length(m) == length(syms) || throw(ArgumentError("got $(length(m)) magnetic moments for $(length(syms)) atoms"))
        new(lattice, syms, Vec3.(frac), Tuple(pbc), String(title), fx, cc, q, m)
    end
end

"""
    with_positions(c, frac; lattice=c.lattice, pbc=c.pbc) -> Crystal

The same atoms (species, constraints, charges, moments, title) at new fractional
positions, optionally in a new cell.
"""
with_positions(c::Crystal, frac; lattice::Lattice=c.lattice, pbc=c.pbc) =
    Crystal(lattice, c.species, frac; pbc, title=c.title, fixed=c.fixed, cell_constraints=c.cell_constraints,
            charges=c.charges, magmoms=c.magmoms)

"""
    with_constraints(c; fixed=c.fixed, cell=c.cell_constraints) -> Crystal

`c` with other atom and/or cell constraints.
"""
with_constraints(c::Crystal; fixed=c.fixed, cell=c.cell_constraints) =
    Crystal(c.lattice, c.species, c.frac; pbc=c.pbc, title=c.title, fixed, cell_constraints=cell, charges=c.charges,
            magmoms=c.magmoms)

"""
    with_charges(c, charges) -> Crystal

`c` with other atomic charges (e; `NaN` or `nothing` = not known).
"""
with_charges(c::Crystal, charges) =
    Crystal(c.lattice, c.species, c.frac; pbc=c.pbc, title=c.title, fixed=c.fixed,
            cell_constraints=c.cell_constraints, charges, magmoms=c.magmoms)

"""Whether any atom has a known charge."""
has_charges(c::Crystal) = any(!isnan, c.charges)

"""
    with_magmoms(c, magmoms) -> Crystal

`c` with other magnetic moments (μB; `NaN` or `nothing` = not known).
"""
with_magmoms(c::Crystal, magmoms) =
    Crystal(c.lattice, c.species, c.frac; pbc=c.pbc, title=c.title, fixed=c.fixed,
            cell_constraints=c.cell_constraints, charges=c.charges, magmoms)

"""Whether any atom has a known magnetic moment."""
has_magmoms(c::Crystal) = any(!isnan, c.magmoms)

"""Whether any atom coordinate or lattice parameter is constrained."""
has_constraints(c::Crystal) = any(any, c.fixed) || c.cell_constraints != FREE_CELL

"""
    crystal_from_cartesian(lattice, species, positions; kwargs...)

Like `Crystal(...)` but takes Cartesian positions (Å) instead of fractional.
"""
function crystal_from_cartesian(lattice::Lattice, species, positions; kwargs...)
    Crystal(lattice, species, [cart_to_frac(lattice, r) for r in positions]; kwargs...)
end

natoms(c::Crystal) = length(c.species)
Base.length(c::Crystal) = natoms(c)

cartesian_positions(c::Crystal) = [frac_to_cart(c.lattice, f) for f in c.frac]

"""Unique species in order of first appearance."""
species_list(c::Crystal) = unique(c.species)

"""Chemical formula like `"Si8"` or `"Ga4As4"` (species in first-appearance order)."""
function formula(c::Crystal)
    join(string(s, count(==(s), c.species)) for s in species_list(c))
end

"""
    wrap(c; atoms=all, tol=1e-8) -> Crystal

Wrap atoms into the cell: map fractional coordinates into [0, 1) along
periodic directions by adding whole lattice vectors (as ASE `wrap`, OVITO
"wrap at periodic boundaries"). The structure itself does not change.
`atoms` limits this to some atom indices; values within `tol` of 1 become 0.
"""
function wrap(c::Crystal; atoms=nothing, tol::Real=1e-8)
    idx = atoms === nothing ? (1:natoms(c)) : atoms
    frac = copy(c.frac)
    for i in idx
        frac[i] = Vec3(ntuple(3) do k
            c.pbc[k] || return frac[i][k]
            x = mod(frac[i][k], 1.0)
            x > 1 - tol ? 0.0 : x
        end)
    end
    with_positions(c, frac)
end

"""
    atoms_outside_cell(c; atoms=all, tol=1e-3) -> Vector{Int}

Atoms with a periodic fractional coordinate outside [0, 1), ignoring atoms
within `tol` of a cell face (those are drawn on the face).
"""
function atoms_outside_cell(c::Crystal; atoms=nothing, tol::Real=1e-3)
    idx = atoms === nothing ? (1:natoms(c)) : atoms
    return [i for i in idx if any(c.pbc[k] && !(-tol <= c.frac[i][k] < 1 + tol) for k in 1:3)]
end

function Base.show(io::IO, c::Crystal)
    print(io, "Crystal(", formula(c), ", ", natoms(c), " atoms)")
end

function Base.show(io::IO, ::MIME"text/plain", c::Crystal)
    println(io, "Crystal ", formula(c), isempty(c.title) ? "" : "  \"$(c.title)\"")
    println(io, "  ", c.lattice)
    @printf(io, "  volume = %.4f Å³, pbc = %s\n", cell_volume(c.lattice), c.pbc)
    n = natoms(c)
    for i in 1:min(n, 10)
        f = c.frac[i]
        @printf(io, "  %-3s %10.6f %10.6f %10.6f\n", c.species[i], f[1], f[2], f[3])
    end
    n > 10 && print(io, "  … ", n - 10, " more atoms")
end
