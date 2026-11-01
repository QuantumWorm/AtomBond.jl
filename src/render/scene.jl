# Put a crystal into a 3D scene.

"""
    draw_structure!(lscene, crystal, settings) -> (DisplayModel, plots, atom_plot)

Draw atoms, bonds and (optionally) the unit cell into an `LScene`. For a
fully non-periodic structure (molecule) the box is drawn as dashed grey lines. Returns the plots so a caller can remove them and redraw without
touching the camera, and the atom plot separately for picking.
"""
function draw_structure!(lscene, c::Crystal, vs::ViewSettings)
    dm = build_display(c, vs)
    plots = Any[]
    if vs.show_cell
        # A molecule's box is drawn dashed grey: exported as the cell, but not periodic.
        append!(plots, any(c.pbc) ? draw_cell!(lscene, c.lattice, vs) : draw_molecule_box!(lscene, c.lattice, vs))
    end
    push!(plots, draw_bonds!(lscene, dm, vs))
    atom_plot = draw_atoms!(lscene, dm, vs)
    push!(plots, atom_plot)
    vs.show_labels && push!(plots, draw_labels!(lscene, dm, vs))
    vs.bond_lengths === :all && vs.show_coordinates && push!(plots, draw_coordinates!(lscene, dm, vs, drawn_atoms(dm, vs)))
    vs.bond_lengths === :all && vs.show_charges && push!(plots, draw_charges!(lscene, dm, vs, drawn_atoms(dm, vs)))
    vs.bond_lengths === :all && vs.show_magmoms && push!(plots, draw_magmoms!(lscene, dm, vs, drawn_atoms(dm, vs)))
    if vs.bond_lengths === :all && bonds_visible(vs)
        vs.show_bond_lengths && push!(plots, draw_bond_lengths!(lscene, dm, vs))
        vs.show_bond_angles && push!(plots, draw_bond_angles!(lscene, dm, vs))
        vs.show_dihedrals && push!(plots, draw_dihedrals!(lscene, dm, vs))
    end
    vs.show_constraints && any(any, c.fixed) && push!(plots, draw_fixed!(lscene, dm, vs))
    filter!(!isnothing, plots)
    return dm, plots, atom_plot
end
