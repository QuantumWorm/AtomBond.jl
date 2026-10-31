# Display settings shared by all render functions.

"""
Display styles, as in Materials Studio:

- `:ball_and_stick` — spheres scaled from covalent radii, cylinder bonds
- `:cpk`            — space-filling spheres at van der Waals radii, no bonds
- `:stick`          — bonds only; atoms are joints the size of the bonds
- `:line`           — bonds as thin lines; only unbonded atoms get a dot
"""
const DISPLAY_STYLES = (:ball_and_stick, :cpk, :stick, :line)

const STYLE_LABELS = Dict(:ball_and_stick => "Ball and stick", :cpk => "CPK (space filling)",
                          :stick => "Stick", :line => "Line")

"""
    ViewSettings(; kwargs...)

Display settings of a viewer (`v.settings`), set by keyword, e.g.
`ViewSettings(style=:cpk, show_bonds=false)`. Change a field of an open
viewer and call [`redraw!`](@ref). Main fields: `style` (see
`DISPLAY_STYLES`), `ball_scale`, `bond_radius`, `show_bonds`, `bond_stubs`,
`show_cell`, `show_boundary_atoms`, `wrap_atoms`, `show_labels`,
`show_triad`, `orthographic`, `bond_lengths` (`:selected` or `:all`),
`show_bond_angles`, `show_dihedrals`, `show_coordinates`, `show_charges`,
`show_magmoms` and `bond_cutoffs` (per element pair).
"""
Base.@kwdef mutable struct ViewSettings
    style::Symbol = :ball_and_stick

    # Atoms
    ball_scale::Float64 = 0.4      # ball & stick: radius = max(ball_scale * r_cov, min_radius)
    min_radius::Float64 = 0.2
    cpk_scale::Float64 = 1.0       # CPK: radius = cpk_scale * r_vdW
    sphere_resolution::Int = 48
    specular::Float64 = 0.4
    shininess::Float64 = 32.0

    # Fixed atoms (constraints for geometry optimization) get a marker:
    # × fully fixed, + some coordinates fixed.
    show_constraints::Bool = true
    # Atoms outside the cell are drawn where they are, unless wrap_atoms.
    wrap_atoms::Bool = false
    # Atoms on cell faces/edges/corners are repeated on the opposite side.
    show_boundary_atoms::Bool = true
    boundary_tol::Float64 = 1e-3   # fractional distance from a face that counts as "on" it

    # Bonds: d ≤ bond_tolerance * (r_cov,i + r_cov,j). Each half is colored by
    # its atom unless bond_color is set.
    show_bonds::Bool = true
    bond_tolerance::Float64 = 1.15
    bond_cutoffs::Dict{Tuple{Symbol,Symbol},Float64} = Dict{Tuple{Symbol,Symbol},Float64}()   # per-pair overrides (kept bonds after substitutions)
    bond_radius::Float64 = 0.12
    bond_color::Union{Nothing,RGBf} = nothing
    bond_stubs::Bool = true        # half-bonds toward atoms outside the cell
    bond_resolution::Int = 24
    line_width::Float64 = 2.5      # :line style

    # Unit cell box; the a, b, c edges from the origin are drawn red/green/blue.
    show_cell::Bool = true
    cell_color::RGBf = RGBf(0.35, 0.35, 0.35)
    cell_linewidth::Float64 = 1.5
    molecule_box_color::RGBf = RGBf(0.6, 0.6, 0.6)   # dashed box of a non-periodic structure
    axis_colors::NTuple{3,RGBf} = (RGBf(0.85, 0.1, 0.1), RGBf(0.1, 0.65, 0.1), RGBf(0.1, 0.25, 0.9))
    show_axis_labels::Bool = false # a/b/c text at the cell edges (the triad is clearer)
    show_triad::Bool = true        # small a/b/c axes in the lower-left corner

    # Labels and selection
    show_labels::Bool = false      # "Si5"-style label on every atom
    label_size::Float64 = 13
    # Bond lengths (Å) at bond midpoints: :selected = bonds between selected
    # atoms (3 or more selected; 2 show the measured distance), :all = every bond.
    bond_lengths::Symbol = :selected
    show_bond_lengths::Bool = true     # the length labels themselves (angles etc. have their own switches)
    bond_length_color::RGBf = RGBf(0.1, 0.25, 0.55)
    # Bond angles (°) inside each angle at its central atom, for the same bonds
    # as the lengths (with :selected: among 4 or more selected atoms).
    show_bond_angles::Bool = false
    bond_angle_color::RGBf = RGBf(0.1, 0.45, 0.15)
    # Dihedral angles (°) i–j–k–l about each bond j–k, on a ring around the bond,
    # for the same bonds (with :selected: among 5 or more selected atoms).
    show_dihedrals::Bool = false
    dihedral_color::RGBf = RGBf(0.45, 0.15, 0.5)
    # Atom coordinates beside the atoms (selected ones, or all with
    # bond_lengths = :all): :fractional, :cartesian (Å), or :auto = fractional
    # for periodic structures, Cartesian for molecules.
    show_coordinates::Bool = false
    coordinates::Symbol = :auto
    coordinate_color::RGBf = RGBf(0.35, 0.25, 0.1)
    # Atomic charges (e) above-left of the atoms, same atoms as the coordinates;
    # only atoms whose charge is known (from the file) get one.
    show_charges::Bool = false
    positive_charge_color::RGBf = RGBf(0.75, 0.1, 0.1)
    negative_charge_color::RGBf = RGBf(0.1, 0.25, 0.75)
    # Magnetic moments (μB) below-left of the atoms, same atoms; only known ones.
    show_magmoms::Bool = false
    spin_up_color::RGBf = RGBf(0.8, 0.4, 0.0)
    spin_down_color::RGBf = RGBf(0.0, 0.5, 0.5)
    highlight_color::RGBAf = RGBAf(1.0, 0.8, 0.0, 0.35)
    measure_color::RGBf = RGBf(0.1, 0.1, 0.1)

    # Scene
    background::RGBf = RGBf(1, 1, 1)
    show_legend::Bool = true
    # Ambient occlusion (darkens crevices). Needs a real GPU: under WSL's
    # software OpenGL (llvmpipe) it just dims the whole image, so it is off by
    # default. It is a window setting, so it must be chosen before opening.
    ssao::Bool = false
    orthographic::Bool = false
    fov::Float64 = 30.0            # perspective field of view (degrees)
    view_direction::Vec3 = Vec3(1.0, -1.2, 0.6)   # from the cell center toward the eye
end

function atom_radius(s::Symbol, vs::ViewSettings)
    vs.style === :cpk && return vs.cpk_scale * vdw_radius(s)
    vs.style === :stick && return vs.bond_radius
    vs.style === :line && return 0.15   # dot for unbonded atoms
    return max(vs.ball_scale * element(s).covalent_radius, vs.min_radius)
end

bonds_visible(vs::ViewSettings) = vs.show_bonds && vs.style !== :cpk
