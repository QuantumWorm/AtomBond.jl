# Conversion to and from AtomsBase systems, the common structure interface of
# the Julia atomistic ecosystem (DFTK, AtomsIO, InteratomicPotentials, ...).
#
#   using DFTK
#   system = to_atomsbase(slab(si, (0, 0, 1)))
#   model  = model_DFT(system; functionals=PBE(), pseudopotentials=...)

"""
    to_atomsbase(c) -> AtomsBase.FlexibleSystem

Cell vectors and Cartesian positions in Å; periodicity from `c.pbc`.
"""
function to_atomsbase(c::Crystal)
    box = Tuple(Vector(v) * u"Å" for v in cell_vectors(c.lattice))
    atoms = [AtomsBase.Atom(s, Vector(r) * u"Å") for (s, r) in zip(c.species, cartesian_positions(c))]
    return AtomsBase.atomic_system(atoms, box, c.pbc)
end

"""
    from_atomsbase(system; title="") -> Crystal

Any AtomsBase system with three cell vectors (e.g. one returned by DFTK or
AtomsIO).
"""
function from_atomsbase(sys; title::AbstractString="")
    cv = AtomsBase.cell_vectors(sys)
    length(cv) == 3 || throw(ArgumentError("need a 3D system with three cell vectors"))
    L = Lattice([ustrip.(u"Å", v) for v in cv]...)
    species = [AtomsBase.atomic_symbol(sys, i) for i in 1:length(sys)]
    positions = [Vec3(ustrip.(u"Å", AtomsBase.position(sys, i))) for i in 1:length(sys)]
    return crystal_from_cartesian(L, species, positions; pbc=Tuple(AtomsBase.periodicity(sys)), title)
end
