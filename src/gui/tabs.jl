# Tabs: several structures in one window.
#
# Each tab is a Document: its structure, the "Original" one, undo/redo
# history, selection, file, 3D view and last oxygen bridge. The active tab's
# state lives in the Viewer's own fields (so all editing code works on it
# unchanged); the other tabs are kept in `v.tabs`, and switching swaps them.
#
# Open…, recent files, New structure and the structure library open a new
# tab, except in place of an untouched empty cell (as `visualize()` starts).
# The row of tabs is above the structure's title: click a tab to show it, ×
# closes it (twice with unsaved changes), + starts a new structure.

const History = Vector{Tuple{Crystal,String,Bool,Union{Nothing,String}}}

mutable struct Document
    crystal::Crystal
    original::Crystal
    selection::Vector{Tuple{Int,Shift}}
    history::History
    future::History
    path::Union{Nothing,String}
    camera::Any        # (eye, lookat, up) of the 3D view, or nothing: frame the structure
    bridge::Any        # the last oxygen bridge (for "Flip last bridge"), or nothing
end

Document(c::Crystal; path=nothing) =
    Document(c, c, Tuple{Int,Shift}[], History(), History(), path === nothing ? nothing : abspath(path), nothing, nothing)

"""A new, empty 5 Å cubic cell (for a new tab)."""
empty_cell() = prototype_crystal("Empty cell", String[], (5.0, 5.0, 5.0))

const TAB_ACTIVE = RGBf(0.80, 0.87, 0.98)
const TAB_IDLE = RGBf(0.92, 0.92, 0.92)
const TAB_MAXNAME = 24

# Tools that work on the structure: they are rebuilt for the new tab when
# tabs switch (New structure and the library do not depend on it).
const STRUCTURE_TOOLS = Dict(:surface => surface_dialog, :export => export_dialog, :atoms => atoms_dialog,
                             :supercell => supercell_dialog, :shift => shift_dialog, :rotate => rotate_dialog,
                             :cell => cell_dialog, :bridge_tool => bridge_dialog)

"""Number of open tabs."""
ntabs(v::Viewer) = length(v.tabs)
"""Number of the tab shown."""
active_tab(v::Viewer) = v.active

"""Save the active tab's state (kept in the viewer) into its Document."""
function stash!(v::Viewer)
    cam = camera(v.lscene)
    v.tabs[v.active] = Document(v.crystal, v.original, copy(v.selection), v.history, v.future, v.path,
                                (cam.eyeposition[], cam.lookat[], cam.upvector[]), get(v.dialogs, :bridge, nothing))
end

"""Show tab `k`: load its Document into the viewer and redraw."""
function show_tab!(v::Viewer, k::Integer)
    d = v.tabs[k]
    v.active = k
    v.crystal, v.original, v.path = d.crystal, d.original, d.path
    v.history, v.future = d.history, d.future
    v.selection = copy(d.selection)
    d.bridge === nothing ? delete!(v.dialogs, :bridge) : (v.dialogs[:bridge] = d.bridge)
    delete!(v.dialogs, :exit_armed)
    redraw!(v)
    v.lattice[] = d.crystal.lattice
    v.header[] = header_text(d.crystal)
    update_legend!(v)
    if d.camera === nothing
        reset_camera!(v)
    else
        eye, lookat, up = d.camera
        Makie.update_cam!(v.lscene.scene, camera(v.lscene), eye, lookat, up)
    end
    update_document!(v)
    build_tabbar!(v)
    # A tool on the right now works on this tab: rebuild it for this structure.
    key = open_tool(v)
    if haskey(STRUCTURE_TOOLS, key)
        close_tool!(v; restore=false)
        STRUCTURE_TOOLS[key](v)
    end
    return v
end

"""
    switch_tab!(viewer, k)

Show the structure of tab `k` (its view, selection and undo history come back).
"""
function switch_tab!(v::Viewer, k::Integer)
    1 <= k <= ntabs(v) || throw(ArgumentError("there is no tab $k"))
    k == v.active && return v
    stash!(v)
    return show_tab!(v, k)
end

"""
    new_tab!(viewer, crystal=empty cell; path=nothing)

Open `crystal` in a new tab and show it. Without a crystal: an empty cell,
with the New structure tool open on the right.
"""
function new_tab!(v::Viewer, c::Crystal=empty_cell(); path=nothing)
    stash!(v)
    push!(v.tabs, Document(c; path))
    show_tab!(v, ntabs(v))
    natoms(c) == 0 && path === nothing && open_tool(v) !== :new && new_structure_dialog(v)
    return v
end

"""An untouched empty cell (as `visualize()` starts): the next structure replaces it instead of opening a tab."""
is_blank(v::Viewer) = v.path === nothing && isempty(v.history) && natoms(v.crystal) == 0

"""
    open_in_tab!(viewer, crystal; path=nothing) -> Bool

Show a structure that was opened, built or taken from the library: in a new
tab, or in place of an untouched empty cell. A file that is already open in
a tab is shown there instead. `path` is the file it came from.
"""
function open_in_tab!(v::Viewer, c::Crystal; path=nothing)
    path === nothing || (path = abspath(path))
    if path !== nothing
        k = findfirst(k -> tab_path(v, k) == path, eachindex(v.tabs))
        if k !== nothing && !tab_modified(v, k)
            switch_tab!(v, k)
            v.status[] = "$(basename(path)) is already open in this tab"
            return true
        end
    end
    if is_blank(v)
        v.tabs[v.active] = Document(c; path)
        show_tab!(v, v.active)
    else
        new_tab!(v, c; path)
    end
    if path !== nothing
        v.saved[path] = c
        add_recent!(path)
        update_document!(v)
    end
    return true
end

# The active tab's state is in the viewer, the others' in their Documents.
tab_crystal(v::Viewer, k) = k == v.active ? v.crystal : v.tabs[k].crystal
tab_path(v::Viewer, k) = k == v.active ? v.path : v.tabs[k].path
tab_history(v::Viewer, k) = k == v.active ? v.history : v.tabs[k].history

"""Unsaved changes in tab `k` (see `is_modified`)."""
function tab_modified(v::Viewer, k)
    path = tab_path(v, k)
    return path === nothing ? !isempty(tab_history(v, k)) : get(v.saved, path, nothing) !== tab_crystal(v, k)
end

"""Name on a tab: the file, else the title or formula; "*" marks unsaved changes."""
function tab_label(v::Viewer, k)
    c, path = tab_crystal(v, k), tab_path(v, k)
    name = path !== nothing ? basename(path) : natoms(c) == 0 ? "Untitled" : isempty(c.title) ? formula(c) : c.title
    length(name) > TAB_MAXNAME && (name = first(name, TAB_MAXNAME - 1) * "…")
    return tab_modified(v, k) ? name * " *" : name
end

"""
    close_tab!(viewer, k=active tab; force=false) -> Bool

Close tab `k`. With unsaved changes the first call only warns in the status
line; a second one (with no edits in between) closes it. Closing the last
tab leaves an empty cell. Returns whether the tab closed.
"""
function close_tab!(v::Viewer, k::Integer=v.active; force::Bool=false)
    1 <= k <= ntabs(v) || throw(ArgumentError("there is no tab $k"))
    c = tab_crystal(v, k)
    if !force && tab_modified(v, k) && get(v.dialogs, :close_armed, nothing) !== c
        v.dialogs[:close_armed] = c
        v.status[] = "Unsaved changes in $(tab_label(v, k)): close the tab again to discard them, or Save (Ctrl+S) first"
        return false
    end
    delete!(v.dialogs, :close_armed)
    name = tab_label(v, k)
    stash!(v)
    deleteat!(v.tabs, k)
    isempty(v.tabs) && push!(v.tabs, Document(empty_cell()))
    if k == v.active
        show_tab!(v, min(k, ntabs(v)))
    else
        k < v.active && (v.active -= 1)
        build_tabbar!(v)
    end
    v.status[] = "Closed $(rstrip(name, ['*', ' ']))"
    return true
end

"""Rebuild the row of tabs (after one is opened, closed or switched to)."""
function build_tabbar!(v::Viewer)
    slot = get(v.dialogs, :tabslot, nothing)
    slot === nothing && return nothing
    old = get(v.dialogs, :tabbar, nothing)
    old === nothing || delete_layout!(old.layout)
    g = Makie.GridLayout(slot[1, 1]; halign=:left, tellwidth=false, colgap=0)
    tabbutton(pos, label, color) = with_click_feedback!(Makie.Button(pos; label, fontsize=13, buttoncolor=color,
                                                                    cornerradius=0, padding=(8, 8, 4, 4)))
    buttons = map(eachindex(v.tabs)) do k
        color = k == v.active ? TAB_ACTIVE : TAB_IDLE
        b = tabbutton(g[1, 3k-2], tab_label(v, k), color)
        Makie.on(_ -> switch_tab!(v, k), b.clicks)
        x = tabbutton(g[1, 3k-1], "×", color)
        Makie.on(_ -> close_tab!(v, k), x.clicks)
        k < ntabs(v) && Makie.Box(g[1, 3k]; width=4, visible=false)       # gap between tabs
        b
    end
    plus = tabbutton(g[1, 3ntabs(v)], "+", TAB_IDLE)
    Makie.on(_ -> new_tab!(v), plus.clicks)
    # An empty last column takes the rest of the row, so the tabs keep their own widths.
    Makie.Box(g[1, 3ntabs(v)+1]; visible=false)
    Makie.colsize!(g, 3ntabs(v) + 1, Makie.Auto(false))
    Makie.colgap!(g, 0)                 # columns added after creation get the default gap
    v.dialogs[:tabbar] = (; layout=g, buttons)
    return g
end

"""Update the names on the tabs (after an edit or a save)."""
function refresh_tab_labels!(v::Viewer)
    bar = get(v.dialogs, :tabbar, nothing)
    bar === nothing && return nothing
    length(bar.buttons) == ntabs(v) || return build_tabbar!(v)
    for (k, b) in enumerate(bar.buttons)
        label = tab_label(v, k)
        b.label[] == label || (b.label[] = label)
    end
end
