# Precompile workload: run what a first session typically does, so the
# compiled code is cached when the package is precompiled and the first
# window and tools open quickly. No OpenGL window is created (`show=false`
# builds the figure only), so this also works without a display.

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    examples = normpath(joinpath(@__DIR__, "..", "examples", "structures"))
    files = filter(f -> structure_format(f) !== nothing,
                   [readdir(examples; join=true); readdir(joinpath(examples, "formats"); join=true)])
    @compile_workload begin
        # Readers, analysis and writers
        for f in files
            try
                c = load_structure(f)
                find_bonds(c)
            catch
            end
        end
        si = load_library("Silicon")
        gan = load_library("GaN")
        water = load_library("Water")
        spacegroup(gan); kpoint_grid(gan)
        sc = supercell(si, 2, 2, 2)
        sl = slab(si, (1, 1, 1); layers=2, vacuum=10)
        fixed = with_constraints(fix_atoms(sl, 1:2); cell=(1, 1, 1, 0, 0, 0))
        for w in (write_qe, write_abinit, write_siesta, write_castep, write_cif, write_xyz, write_poscar)
            sprint(io -> w(io, fixed))
        end
        nanotube(6, 5)

        # The window (figure only), its panel and tools, and common actions
        v = visualize(gan; show=false)
        for tool in (library_dialog, supercell_dialog, surface_dialog, export_dialog, cell_dialog, atoms_dialog,
                     shift_dialog, rotate_dialog, bridge_dialog, nanotube_dialog, new_structure_dialog)
            tool(v)
        end
        close_tool!(v)
        select!(v, 1); select!(v, 2; add=true); clear_selection!(v)
        edit!(c -> supercell(c, 2, 2, 1), v, "make supercell")
        undo!(v)
        view_along!(v, :c); rotate_view!(v, 90); reset_camera!(v)
        set_bond_lengths!(v, :all)
        open_in_tab!(v, water)
        set_crystal!(v, sc)
        for name in ("File", "Edit", "Build", "View")
            menubar_items(v, name)
        end
    end
end
