# Atoms as instanced spheres (one meshscatter call for the whole crystal).

"""
The atoms actually drawn as spheres, in plot order (so a pick index on the
atom plot maps straight into this list). Line style draws only unbonded atoms.
"""
function drawn_atoms(dm::DisplayModel, vs::ViewSettings)
    vs.style === :line || return dm.atoms
    bonded = Set(i for b in dm.bonds for i in (b.i, b.j))
    return filter(a -> a.index ∉ bonded, dm.atoms)
end

function draw_atoms!(ax, dm::DisplayModel, vs::ViewSettings)
    atoms = drawn_atoms(dm, vs)
    isempty(atoms) && return nothing
    species = [dm.crystal.species[a.index] for a in atoms]
    positions = [Point3f(a.position) for a in atoms]
    radii = Float32[atom_radius(s, vs) for s in species]
    colors = [element(s).color for s in species]
    sphere = Makie.Tessellation(Sphere(Point3f(0), 1.0f0), mesh_detail(vs, length(dm.atoms))[1])
    Makie.meshscatter!(ax, positions;
        marker=sphere, markersize=radii, color=colors,
        specular=vs.specular, shininess=vs.shininess, ssao=vs.ssao)
end
