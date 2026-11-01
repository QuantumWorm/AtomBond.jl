# Main window.

"""
Everything a window needs to redraw itself when a setting changes.
"""
mutable struct Viewer
    crystal::Crystal
    original::Crystal                       # as loaded, for "Original" in the panel
    settings::ViewSettings
    figure::Makie.Figure
    lscene::Makie.LScene
    display::DisplayModel
    plots::Vector{Any}
    atom_plot::Any                          # sphere plot, for picking
    selection::Vector{Tuple{Int,Shift}}     # (atom index, image) in click order
    overlay::Vector{Any}                    # selection highlight + measurement plots
    info::Makie.Observable{String}          # description of the selection
    triad::Union{Nothing,Makie.Scene}
    lattice::Makie.Observable{Lattice}      # drives the triad
    header::Makie.Observable{String}
    status::Makie.Observable{String}
    screen::Union{Nothing,GLMakie.Screen}   # the open window, if any
    legend::Any
    history::Vector{Tuple{Crystal,String,Bool,Union{Nothing,String}}}   # undo stack: (crystal, what, reframe, file)
    future::Vector{Tuple{Crystal,String,Bool,Union{Nothing,String}}}    # redo stack
    textboxes::Vector{Any}                        # panel text boxes (shortcuts pause while typing)
    dialogs::Dict{Symbol,Any}                     # GUI parts and state by name (:tool on the right, :panel, :menu, :tabbar, …)
    folder::String                                # where the Open tool starts
    path::Union{Nothing,String}                   # file of the structure (Save writes here), if any
    saved::Dict{String,Crystal}                   # file => structure as last opened/saved there, for "modified"
    document::Makie.Observable{String}            # file name + "modified", in the menu bar
    tabs::Vector{Any}                             # one Document per tab (tabs.jl); the active one's state is above
    active::Int                                   # the tab shown
end

"""Replace the structure plots after a settings change; the camera stays put."""
function redraw!(v::Viewer)
    foreach(p -> delete!(v.lscene.scene, p), v.plots)
    v.display, v.plots, v.atom_plot = draw_structure!(v.lscene, v.crystal, v.settings)
    update_selection!(v)    # positions/radii may have changed
    return v
end

"""
    set_crystal!(viewer, crystal; reframe=true, keep_selection=false)

Show a different structure in the same window. This is not recorded for
undo; edits from the panel go through [`edit!`](@ref). The selection is
cleared unless `keep_selection` (for edits that keep atom numbering).
"""
function set_crystal!(v::Viewer, c::Crystal; reframe::Bool=true, keep_selection::Bool=false)
    species_changed = species_list(c) != species_list(v.crystal)
    v.crystal = c
    keep_selection ? filter!(k -> k[1] <= natoms(c), v.selection) : empty!(v.selection)
    redraw!(v)
    v.lattice[] = c.lattice
    v.header[] = header_text(c)
    update_document!(v)
    species_changed && update_legend!(v)
    reframe && reset_camera!(v)
    return v
end

# Status messages that ask for input or report a failure are shown in red.
const PROBLEM_PREFIXES = ("Type ", "Unsaved ", "Enter ", "Could not", "Cannot", "Select ", "Unknown", "Nothing to",
                          "Choose ", "[0 0 0]")
status_color(s) = any(p -> startswith(s, p), PROBLEM_PREFIXES) ? RGBf(0.75, 0.1, 0.1) : RGBf(0.3, 0.3, 0.3)

"""Element legend under the 3D view; rebuilt when the set of elements changes."""
function update_legend!(v::Viewer)
    v.legend === nothing || delete!(v.legend)
    v.legend = nothing
    v.settings.show_legend || return
    spc = species_list(v.crystal)
    isempty(spc) && return
    markers = [Makie.MarkerElement(marker=:circle, color=element(s).color,
                                   strokecolor=:black, strokewidth=0.5, markersize=16) for s in spc]
    v.legend = Makie.Legend(v.figure[3, 1], markers, string.(spc); orientation=:horizontal,
                            framevisible=false, tellheight=true, tellwidth=false)
end

"""Title line: formula, space group (when periodic), and the file's title."""
function header_text(c::Crystal)
    parts = [natoms(c) == 0 ? "No atoms" : formula(c)]
    sg = try
        spacegroup(c)
    catch
        nothing
    end
    sg === nothing || push!(parts, "$(sg[2]) ($(sg[1]))")
    isempty(c.title) || push!(parts, c.title)
    return join(parts, "  —  ")
end

view_along!(v::Viewer, which) = view_along!(v.lscene, v.display, v.settings, which)
rotate_view!(v::Viewer, angle::Real) = (rotate_view!(v.lscene, angle); v)
reset_camera!(v::Viewer) = reset_camera!(v.lscene, v.display, v.settings)

"""
    visualize(crystal_or_path; settings=ViewSettings(), size=(1250, 960), panel=true, show=true) -> Viewer

Open a GLMakie window showing the structure, with a control panel on the
right (`panel=false` hides it). Rotate with left-drag, zoom with the scroll
wheel, pan with right-drag.
"""
function visualize(c::Crystal; settings::ViewSettings=ViewSettings(),
                   size=(1250, 990), panel::Bool=true, show::Bool=true, folder::AbstractString=pwd(),
                   path::Union{Nothing,AbstractString}=nothing)
    vs = settings
    GLMakie.activate!()
    fig = Makie.Figure(; size, backgroundcolor=vs.background)

    header = Makie.Observable(header_text(c))
    top = Makie.GridLayout(fig[1, 1]; tellwidth=false, rowgap=6)     # row of tabs (with the panel), title
    Makie.Label(top[panel ? 2 : 1, 1], header; fontsize=16, halign=:left, tellwidth=false)

    # SSAO radius/bias are in scene units (Å): sized for atoms ~0.2–2 Å.
    lscene = Makie.LScene(fig[2, 1]; show_axis=false,
        scenekw=(ssao=Makie.SSAO(radius=1.2, bias=0.05, blur=2),))
    setup_camera!(lscene, vs)
    dm, plots, atom_plot = draw_structure!(lscene, c, vs)
    reset_camera!(lscene, dm, vs)
    lattice = Makie.Observable(c.lattice)
    # The a b c axes indicator, also for molecules (the axes of their box).
    triad = vs.show_triad ? add_triad!(fig.scene, lscene, lattice, vs) : nothing

    v = Viewer(c, c, vs, fig, lscene, dm, plots, atom_plot, Tuple{Int,Shift}[], Any[],
               Makie.Observable(SELECTION_HINT), triad, lattice, header, Makie.Observable(""), nothing,
               nothing, [], [], [], Dict{Symbol,Any}(), folder, path === nothing ? nothing : abspath(path),
               path === nothing ? Dict{String,Crystal}() : Dict(abspath(path) => c), Makie.Observable(""),
               Any[Document(c; path)], 1)
    update_document!(v)
    update_legend!(v)
    install_picking!(v)
    install_ripple!(v)          # before the popup menu's scene, so menus draw on top
    install_context_menu!(v)
    install_shortcuts!(v)
    if panel
        build_menubar!(fig[0, 1:2], v)
        v.dialogs[:has_panel] = true
        show_panel!(v)
        v.dialogs[:tabslot] = top
        build_tabbar!(v)
        # Footer: what is selected / measured, then messages from the panel.
        Makie.Label(fig[4, 1:2], v.info; halign=:left, fontsize=14, tellwidth=false)
        Makie.Label(fig[5, 1:2], v.status; halign=:left, fontsize=13, color=Makie.lift(status_color, v.status),
                    tellwidth=false)
    end

    if show
        # SSAO is a per-screen option that GLMakie's save/colorbuffer would not
        # carry over, so the viewer always renders through a screen of its own.
        v.screen = GLMakie.Screen(; title="AtomBond — $(formula(c))", ssao=vs.ssao)
        display(v.screen, fig)
    end
    return v
end

visualize(path::AbstractString; kwargs...) =
    visualize(load_structure(path); folder=dirname(abspath(path)), path, kwargs...)

"""
    visualize(; kwargs...) -> Viewer

Start without a structure: an empty 5 Å cubic cell, with the New structure
window open (a prototype such as diamond or wurtzite with your elements, or
an empty cell to fill with Add). File → Open… and the library work as usual.
"""
function visualize(; show::Bool=true, kwargs...)
    v = visualize(prototype_crystal("Empty cell", String[], (5.0, 5.0, 5.0)); show, kwargs...)
    show && new_structure_dialog(v)
    return v
end

"""
Render the whole figure to an image, using the open window or an offscreen
screen. Images get full sphere and bond detail even when the interactive view
uses less (`ViewSettings.adaptive_detail`).
"""
function render_image(v::Viewer)
    vs = v.settings
    reduced = mesh_detail(vs, length(v.display.atoms)) != (vs.sphere_resolution, vs.bond_resolution)
    reduced || return render_image_now(v)
    vs.adaptive_detail = false
    try
        redraw!(v)
        return render_image_now(v)
    finally
        vs.adaptive_detail = true
        redraw!(v)
    end
end

function render_image_now(v::Viewer)
    v.screen === nothing || return Makie.colorbuffer(v.screen)
    # Already drawn offscreen (a figure can be on one GLMakie screen only).
    s = Makie.getscreen(v.figure.scene)
    s isa GLMakie.Screen && return Makie.colorbuffer(s)
    screen = GLMakie.Screen(; visible=false, ssao=v.settings.ssao)
    try
        display(screen, v.figure)
        return Makie.colorbuffer(screen)
    finally
        close(screen)
    end
end

"""
    save_view(path, viewer)

Save just the 3D view (no panel, header or legend) as an image.
"""
function save_view(path::AbstractString, v::Viewer)
    img = render_image(v)
    # Crop the figure image to the 3D view. The viewport is in figure units
    # with the origin at the bottom left; the image may be scaled by px_per_unit.
    s = size(img, 2) / Makie.widths(v.figure.scene.viewport[])[1]
    vp = v.lscene.scene.viewport[]
    (x, y), (w, h) = round.(Int, s .* Makie.origin(vp)), round.(Int, s .* Makie.widths(vp))
    H = size(img, 1)
    Makie.save(path, img[H-(y+h)+1:H-y, x+1:x+w])
    return path
end

"""
    snapshot(outfile, crystal_or_path; kwargs...)

Render to an image file (PNG) without opening a window.
"""
function snapshot(outfile::AbstractString, c; kwargs...)
    v = visualize(c; show=false, panel=false, kwargs...)
    Makie.save(outfile, render_image(v))
    return outfile
end
