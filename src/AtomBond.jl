"""
AtomBond — build, edit and view crystal structures (ball-and-stick, GLMakie).

Layers:
- `core/`     Element data, Lattice, Crystal (no graphics)
- `io/`       File formats
- `render/`   Draw a Crystal into a Makie scene
- `gui/`      Windows and interaction
"""
module AtomBond

using LinearAlgebra
import Dates
import TOML
using Printf
using StaticArrays
using GeometryBasics: Point3f, Vec2f, Vec3f, Sphere, Cylinder
using NearestNeighbors: KDTree, inrange
import Spglib
import AtomsBase
using Unitful: @u_str, ustrip
import GLMakie
import GLMakie: Makie
using GLMakie.Makie: RGBf, RGBAf

include("core/elements.jl")
include("core/lattice.jl")
include("core/crystal.jl")

include("io/poscar.jl")
include("io/xyz.jl")
include("io/cif.jl")
include("io/castep.jl")
include("io/xsf.jl")
include("io/qe.jl")
include("io/codes.jl")
include("io/vasp_out.jl")
include("io/molecular.jl")
include("io/msi.jl")
include("io/export.jl")
include("io/atomsbase.jl")
include("io/load.jl")
include("io/library.jl")

include("analysis/neighbors.jl")
include("analysis/bonds.jl")
include("analysis/symmetry.jl")
include("analysis/geometry.jl")

include("build/supercell.jl")
include("build/edit.jl")
include("build/surface.jl")
include("build/bridge.jl")
include("build/prototypes.jl")
include("build/nanotube.jl")
include("build/constraints.jl")

include("render/style.jl")
include("render/display.jl")
include("render/atoms.jl")
include("render/bonds.jl")
include("render/overlay.jl")
include("render/cell.jl")
include("render/camera.jl")
include("render/triad.jl")
include("render/scene.jl")

include("gui/app.jl")
include("gui/panel.jl")
include("gui/picking.jl")
include("gui/feedback.jl")
include("gui/contextmenu.jl")
include("gui/menubar.jl")
include("gui/history.jl")
include("gui/dialogs.jl")
include("gui/filebrowser.jl")
include("gui/tabs.jl")

# core
export Element, element, element_property, vdw_radius, ismetal, atomic_mass, override_elements!, use_vesta_colors!, Lattice, Crystal, crystal_from_cartesian, cubic
export cell_vectors, cell_volume, lattice_parameters, frac_to_cart, cart_to_frac
export natoms, cartesian_positions, species_list, formula, wrap, atoms_outside_cell
export with_positions, with_constraints, has_constraints, fix_atoms, free_atoms, periodic_box, as_molecule
export describe_cell_constraints
# io
export read_poscar, write_poscar, read_xyz, write_xyz, read_cif, load_structure, save_structure
export read_castep_cell, read_castep_geom, read_castep_geom_frames, read_castep_output, structure_format
export with_charges, has_charges, with_magmoms, has_magmoms
export read_xsf, write_xsf, read_cube, read_qe, read_qe_output, read_aims, write_aims, read_siesta, read_abinit
export read_xdatcar, read_vasprun, read_pdb, write_pdb, read_lammps_data, write_lammps_data
export read_car, write_car, read_xsd, read_res, write_res
export kpoint_grid, write_kpoints, write_vasp_inputs, write_qe, write_castep, write_cif, export_inputs
export to_atomsbase, from_atomsbase
# analysis
export Bond, find_bonds, coordination, neighbor_pairs, bond_cutoff
export spacegroup, symmetry_dataset, primitive_cell, conventional_cell
export distance, bond_angle, dihedral, measure, min_image_distance
# build
export supercell, delete_atoms, substitute, add_atom, move_atoms, translate_atoms, rotate_atoms, rotation_matrix, closest_contact, set_position, replace_species
export surface_cell, slab, best_termination, add_vacuum, stack_cells, superlattice
export prototype_crystal, nanotube
export bridge_position, add_bridge
# render / gui
export ViewSettings, DISPLAY_STYLES, draw_structure!
export view_direction, view_along!, reset_camera!
export Viewer, visualize, redraw!, set_crystal!, save_view, snapshot
export select!, select_all!, clear_selection!, edit!, undo!, redo!, open_tool, close_tool!, click_feedback!
export set_bond_lengths!, set_bond_angles!, set_dihedrals!, set_coordinates!, set_charge_labels!, set_magmom_labels!, set_bond_length_labels!
export LibraryEntry, structure_library, library_entry, load_library, open_library!
export open_structure!, save!, save_to!, exit!
export new_tab!, switch_tab!, close_tab!, open_in_tab!, ntabs, active_tab

"""
In WSL2, Mesa's GPU OpenGL driver (D3D12) deadlocks when many plots are
added to an open window, as when a tool opens in the right-hand column. Use
software rendering there (no window is open yet at load time), unless
ATOMBOND_GPU=1 or the user already chose a driver with LIBGL_ALWAYS_SOFTWARE
or GALLIUM_DRIVER. Returns whether it switched.
"""
function wsl_software_rendering!(env=ENV)
    Sys.islinux() || return false
    get(env, "ATOMBOND_GPU", "") == "1" && return false
    (haskey(env, "LIBGL_ALWAYS_SOFTWARE") || haskey(env, "GALLIUM_DRIVER")) && return false
    wsl = try
        occursin("microsoft", lowercase(read("/proc/version", String)))
    catch
        false
    end
    wsl || return false
    env["LIBGL_ALWAYS_SOFTWARE"] = "1"
    return true
end

function __init__()
    load_user_elements!()      # user overrides of element colors/radii (~/.atombond/elements.csv)
    wsl_software_rendering!()
end

end # module
