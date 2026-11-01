# Conversion to and from AtomsBase systems, the common structure interface of
# the Julia atomistic ecosystem (DFTK, AtomsIO, InteratomicPotentials, ...).
#
#   using DFTK
#   system = to_atomsbase(slab(si, (0, 0, 1)))
#   model  = model_DFT(system; functionals=PBE(), pseudopotentials=...)

"""
    to_atomsbase(c) -> AtomsBase.FlexibleSystem

Cell vectors and Cartesian positions in Å; periodicity from `c.pbc`. Known
initial magnetic moments become the atoms' `magnetic_moment` (μB, read by
DFTK), and fully fixed atoms get `clamp=true` (kept fixed by
GeometryOptimization.jl).
"""
function to_atomsbase(c::Crystal)
    box = Tuple(Vector(v) * u"Å" for v in cell_vectors(c.lattice))
    atoms = map(enumerate(zip(c.species, cartesian_positions(c)))) do (i, (s, r))
        props = Pair{Symbol,Any}[]
        isnan(c.magmoms[i]) || push!(props, :magnetic_moment => c.magmoms[i])
        all(c.fixed[i]) && push!(props, :clamp => true)
        AtomsBase.Atom(s, Vector(r) * u"Å"; props...)
    end
    return AtomsBase.atomic_system(atoms, box, c.pbc)
end

"""
    from_atomsbase(system; title="") -> Crystal

Any AtomsBase system with three cell vectors (e.g. one returned by DFTK,
GeometryOptimization or AtomsIO). Atoms' `magnetic_moment` (a number, μB)
and `clamp=true` (fully fixed) are kept.
"""
function from_atomsbase(sys; title::AbstractString="")
    cv = AtomsBase.cell_vectors(sys)
    length(cv) == 3 || throw(ArgumentError("need a 3D system with three cell vectors"))
    L = Lattice([ustrip.(u"Å", v) for v in cv]...)
    species = [AtomsBase.atomic_symbol(sys, i) for i in 1:length(sys)]
    positions = [Vec3(ustrip.(u"Å", AtomsBase.position(sys, i))) for i in 1:length(sys)]
    c = crystal_from_cartesian(L, species, positions; pbc=Tuple(AtomsBase.periodicity(sys)), title)
    moments = [atom_property(sys, i, :magnetic_moment) for i in 1:length(sys)]
    any(m -> m isa Real, moments) && (c = with_magmoms(c, [m isa Real ? Float64(m) : NaN for m in moments]))
    clamped = [atom_property(sys, i, :clamp) === true for i in 1:length(sys)]
    any(clamped) && (c = fix_atoms(c, findall(clamped)))
    return c
end

"""Property `key` of atom `i` of an AtomsBase system, or `nothing`."""
function atom_property(sys, i, key)
    atom = sys[i]
    return haskey(atom, key) ? atom[key] : nothing
end
