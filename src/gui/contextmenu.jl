# Popup menus in the main window: the right-click menu for atoms, and the
# drop-down menus of the menu bar (menubar.jl).
#
#   right-click on an atom      select it (unless it is already selected) and
#                               open the menu for the selection
#   right-click on empty space  open the menu if atoms are selected
#   right-drag                  still pans the view (only a click opens the menu)
#
# Makie has no popup menus, so a menu is drawn in a pixel-space child scene
# on top of the figure. While it is open it takes the mouse and Esc: a click
# on an item runs it, a click anywhere else only closes the menu.

const MENU_WIDTH = 220
const MENU_ROW = 24
const MENU_SEPARATOR = 9
const MENU_FONTSIZE = 14
const MENU_PRIORITY = 200       # above picking (10) and the camera

"""
One line of a menu. `kind` is :item (runs `action`), :title (bold grey
heading), :disabled (grey, does nothing) or :separator (a line).
`shortcut` is shown on the right.
"""
struct MenuItem
    label::String
    action::Union{Nothing,Function}
    shortcut::String
    kind::Symbol
end
MenuItem(label, action; shortcut="") = MenuItem(label, action, shortcut, action === nothing ? :title : :item)
"""An item that is shown grey and does nothing when `enabled` is false."""
menu_item(label, action, enabled::Bool=true; shortcut="") =
    MenuItem(label, action, shortcut, enabled ? :item : :disabled)
const SEPARATOR = MenuItem("", nothing, "", :separator)

row_height(it::MenuItem) = it.kind === :separator ? MENU_SEPARATOR : MENU_ROW

# One set of plots draws whichever menu is open; opening a menu only updates
# their data. (Deleting and re-adding plots left stale ones on screen in
# GLMakie, covering the highlight of the next menu.)
mutable struct ContextMenu
    scene::Makie.Scene
    items::Vector{MenuItem}
    hover::Makie.Observable{Int}        # row under the mouse, 0 = none
    isopen::Bool
    owner::Any                          # what opened it (a menu bar title), or nothing
    width::Int
    height::Int
    background::Makie.Observable{Makie.Rect2f}
    lines::Makie.Observable{Vector{Makie.Point2f}}                  # separators + border, as segments
    texts::Makie.Observable{Vector{Tuple{String,Makie.Point2f,Symbol}}}    # labels: (text, position, kind)
    shortcuts::Makie.Observable{Vector{Tuple{String,Makie.Point2f}}}
    flash::Makie.Observable{Bool}       # chosen item shown in the strong colour
    generation::Int                     # counts openings (a delayed hide must not hide a newer menu)
end

"""All atoms of element `el` (home-cell copies) become the selection."""
function select_species!(v::Viewer, el::Symbol)
    v.selection = [(i, Shift(0, 0, 0)) for i in eachindex(v.crystal.species) if v.crystal.species[i] == el]
    update_selection!(v)
end

"""Menu items for the current selection (the right-click menu)."""
function menu_items(v::Viewer)
    c = v.crystal
    idx = selected_indices(v)
    isempty(idx) && return MenuItem[]
    spc = unique(c.species[idx])
    title = length(idx) == 1 ? display_label(c, only(selected_atoms(v))) * "  ($(c.species[idx[1]]))" :
            "$(length(idx)) atoms: " * join(("$s×$(count(==(s), c.species[idx]))" for s in spc), ", ")
    items = [MenuItem(title, nothing),
             MenuItem("Edit atoms…", () -> atoms_dialog(v)),
             MenuItem("Shift atoms…", () -> shift_dialog(v)),
             MenuItem("Rotate atoms…", () -> rotate_dialog(v))]
    # Quick element changes to the other species already in the structure.
    for el in species_list(c)
        spc == [el] && continue
        push!(items, MenuItem("Change to $el", () -> substitute_selected!(v, el)))
    end
    if length(selected_atoms(v)) == 3
        push!(items, MenuItem("Add oxygen bridge", () -> add_bridge!(v; bridge_settings(v)...); shortcut="1–2, side 3"))
    end
    last = get(v.dialogs, :bridge, nothing)
    last !== nothing && last.crystal === c && push!(items, MenuItem("Flip oxygen bridge", () -> flip_bridge!(v)))
    nout = length(atoms_outside_cell(c; atoms=idx))
    nout > 0 && push!(items, MenuItem("Wrap $(plural(nout, "atom")) into cell", () -> wrap_atoms!(v; atoms=idx)))
    length(spc) == 1 && push!(items, MenuItem("Select all $(spc[1])", () -> select_species!(v, spc[1])))
    fixed = count(i -> any(c.fixed[i]), idx)
    fixed < length(idx) && push!(items, MenuItem("Fix " * plural(length(idx), "atom"), () -> fix_selected!(v);
                                                 shortcut="constraint"))
    fixed > 0 && push!(items, MenuItem("Free " * plural(fixed, "fixed atom"), () -> free_selected!(v)))
    push!(items, MenuItem("Delete " * plural(length(idx), "atom"), () -> delete_selected!(v); shortcut="Del"),
                 MenuItem("Clear selection", () -> clear_selection!(v); shortcut="Esc"))
    return items
end

function menu_scene(v::Viewer)
    haskey(v.dialogs, :menu) && return v.dialogs[:menu]
    # Created on first use, after everything else, so it is drawn on top.
    scene = Makie.Scene(v.figure.scene; viewport=Makie.Rect2i(0, 0, 1, 1), clear=false, visible=false)
    Makie.campixel!(scene)
    P = Makie.Point2f
    m = ContextMenu(scene, MenuItem[], Makie.Observable(0), false, nothing, MENU_WIDTH, 1,
                    Makie.Observable(Makie.Rect2f(0, 0, 1, 1)), Makie.Observable([P(0, 0), P(0, 0)]),
                    Makie.Observable([("", P(0, 0), :item)]), Makie.Observable([("", P(0, 0))]),
                    Makie.Observable(false), 0)
    # GLMakie draws lines and text in a later pass than meshes, for all scenes,
    # so the menu's background must be in that pass too (fxaa=false), and
    # overdraw ignores the depth of the 3D view.
    top = (fxaa=false, overdraw=true)
    highlight = Makie.lift(m.hover) do k
        (k == 0 || k > length(m.items)) && return Makie.Rect2f(0, 0, 0, 0)
        t, b = row_extents(m.items)[k]
        Makie.Rect2f(1, m.height - b, m.width - 2, b - t)
    end
    Makie.poly!(scene, m.background; color=Makie.RGBf(0.97, 0.97, 0.97), top...)
    Makie.poly!(scene, highlight; color=Makie.lift(f -> f ? FLASH_COLOR : Makie.RGBf(0.80, 0.87, 0.98), m.flash),
                top...)
    Makie.linesegments!(scene, m.lines; color=:gray65, linewidth=1, top...)
    # Titles bold, other labels regular (one text plot each; a text plot has one font).
    for (bold, font) in ((false, :regular), (true, :bold))
        sel = Makie.lift(t -> (u = filter(x -> (x[3] === :title) == bold, t); isempty(u) ? [("", Makie.Point2f(0, 0), :item)] : u),
                         m.texts)
        Makie.text!(scene, Makie.lift(t -> getindex.(t, 2), sel); text=Makie.lift(t -> getindex.(t, 1), sel),
                    color=Makie.lift(t -> [k === :item ? Makie.RGBf(0, 0, 0) : Makie.RGBf(0.5, 0.5, 0.5) for (_, _, k) in t], sel),
                    font, fontsize=MENU_FONTSIZE, align=(:left, :center), top...)
    end
    Makie.text!(scene, Makie.lift(t -> last.(t), m.shortcuts); text=Makie.lift(t -> first.(t), m.shortcuts),
                fontsize=12, align=(:right, :center), color=:gray45, top...)
    v.dialogs[:menu] = m
    return m
end

# Top and bottom (pixels from the menu's top edge) of each row.
function row_extents(items)
    tops = cumsum([0; [row_height(it) for it in items]]) .+ 2
    return [(tops[k], tops[k+1]) for k in eachindex(items)]
end

"""
    open_menu!(viewer, items, xy; width=MENU_WIDTH, owner=nothing) -> ContextMenu or nothing

Show a popup menu with its top-left corner at window pixel `xy` (moved
inside the window if needed). Nothing happens for an empty item list.
"""
function open_menu!(v::Viewer, items::Vector{MenuItem}, xy; width::Integer=MENU_WIDTH, owner=nothing)
    isempty(items) && return nothing
    m = menu_scene(v)
    close_context_menu!(v)
    ext = row_extents(items)
    h = ext[end][2] + 2
    m.items, m.owner, m.width, m.height = items, owner, width, h
    W, H = v.figure.scene.viewport[].widths
    x = clamp(round(Int, xy[1]), 0, max(0, W - width))
    y = clamp(round(Int, xy[2]) - h, 0, max(0, H - h))
    m.scene.viewport[] = Makie.Rect2i(x, y, width, h)

    P = Makie.Point2f
    lines = P[(0.5, 0.5), (width - 0.5, 0.5), (width - 0.5, 0.5), (width - 0.5, h - 0.5),
              (width - 0.5, h - 0.5), (0.5, h - 0.5), (0.5, h - 0.5), (0.5, 0.5)]
    texts, shortcuts = Tuple{String,P,Symbol}[], Tuple{String,P}[]
    for (k, it) in enumerate(items)
        ymid = h - (ext[k][1] + ext[k][2]) / 2
        if it.kind === :separator
            push!(lines, P(6, ymid), P(width - 6, ymid))
            continue
        end
        push!(texts, (it.label, P(10, ymid), it.kind))
        isempty(it.shortcut) || push!(shortcuts, (it.shortcut, P(width - 10, ymid)))
        # A line under a title that starts the menu.
        k == 1 && it.kind === :title && push!(lines, P(6, h - ext[k][2]), P(width - 6, h - ext[k][2]))
    end
    m.background[] = Makie.Rect2f(0, 0, width, h)
    m.lines[] = lines
    m.texts[] = texts
    m.shortcuts[] = isempty(shortcuts) ? [("", P(0, 0))] : shortcuts
    m.hover[] = 0
    m.flash[] = false
    m.generation += 1
    m.isopen = true
    m.scene.visible[] = true
    return m
end

"""
    open_context_menu!(viewer, xy) -> ContextMenu or nothing

Show the right-click menu for the selected atoms near window pixel `xy`.
Nothing happens with no selection.
"""
open_context_menu!(v::Viewer, xy) =
    # Down and to the right, a little off the cursor so the atom stays visible.
    open_menu!(v, menu_items(v), (xy[1] + 12, xy[2] - 12))

"""Close the open popup menu; with `linger` (s) it stays visible that long (inactive)."""
function close_context_menu!(v::Viewer; linger::Real=0)
    m = get(v.dialogs, :menu, nothing)
    (m === nothing || !m.isopen) && return false
    m.isopen = false
    hide() = (m.scene.visible[] = false; m.scene.viewport[] = Makie.Rect2i(0, 0, 1, 1))
    if linger > 0
        g = m.generation
        after_frames(m.scene; seconds=linger) do
            m.generation == g && !m.isopen && hide()
        end
    else
        hide()
    end
    return true
end

menu_isopen(v::Viewer) = (m = get(v.dialogs, :menu, nothing)) !== nothing && m.isopen

"""Menu row at window pixel `xy` (0 = inside but not on an active item, `nothing` = outside)."""
function menu_row_at(v::Viewer, xy)
    m = v.dialogs[:menu]
    vp = m.scene.viewport[]
    o, w = Makie.origin(vp), Makie.widths(vp)
    x, y = xy[1] - o[1], xy[2] - o[2]
    (0 <= x <= w[1] && 0 <= y <= w[2]) || return nothing
    from_top = w[2] - y
    k = findfirst(e -> e[1] <= from_top < e[2], row_extents(m.items))
    return k !== nothing && m.items[k].kind === :item ? k : 0
end

"""Window pixel at the centre of menu row `k` (for tests)."""
function menu_row_pixel(v::Viewer, k)
    m = v.dialogs[:menu]
    vp = m.scene.viewport[]
    o, w = Makie.origin(vp), Makie.widths(vp)
    t, b = row_extents(m.items)[k]
    return (o[1] + w[1] / 2, o[2] + w[2] - (t + b) / 2)
end

"""Run menu row `k` (closing the menu first)."""
function choose_menu_item!(v::Viewer, k)
    m = v.dialogs[:menu]
    it = m.items[k]
    # Flash the chosen item, then let the menu disappear (the action runs now).
    CLICK_FLASH[] && it.kind === :item && (m.hover[] = k; m.flash[] = true)
    close_context_menu!(v; linger=CLICK_FLASH[] ? FLASH_TIME : 0)
    it.kind === :item && (click_sound(); it.action())
    return nothing
end

"""Right-click (press and release without moving) in the 3D view."""
function handle_right_click!(v::Viewer)
    a = atom_under_mouse(v)
    if a !== nothing && (a.index, a.image) ∉ v.selection
        select!(v, a.index, a.image)
    end
    click_ripple!(v, Makie.events(v.figure.scene).mouseposition[], :menu)
    open_context_menu!(v, Makie.events(v.figure.scene).mouseposition[])
end

function install_context_menu!(v::Viewer)
    ev = Makie.events(v.figure.scene)
    Makie.on(ev.mouseposition; priority=MENU_PRIORITY) do xy
        menu_isopen(v) || return Makie.Consume(false)
        m = v.dialogs[:menu]
        # Moving onto another menu bar title switches to its menu.
        t = m.owner === nothing ? nothing : menubar_title_at(v, xy)
        t === nothing || t == m.owner || (open_menubar_menu!(v, t); return Makie.Consume(false))
        k = something(menu_row_at(v, xy), 0)
        m.hover[] == k || (m.hover[] = k)
        return Makie.Consume(false)
    end
    # While open, the menu takes every button press: an item runs, elsewhere closes.
    Makie.on(ev.mousebutton; priority=MENU_PRIORITY) do e
        menu_isopen(v) || return Makie.Consume(false)
        e.action == Makie.Mouse.press || return Makie.Consume(true)
        k = menu_row_at(v, ev.mouseposition[])
        t = k === nothing ? menubar_title_at(v, ev.mouseposition[]) : nothing
        if e.button == Makie.Mouse.left && k !== nothing && k > 0
            choose_menu_item!(v, k)
        elseif t !== nothing && t != v.dialogs[:menu].owner
            open_menubar_menu!(v, t)         # another menu bar title
        elseif k === nothing
            close_context_menu!(v)           # (a click on its own title closes it too)
        end
        return Makie.Consume(true)
    end
    Makie.on(ev.keyboardbutton; priority=MENU_PRIORITY) do e
        menu_isopen(v) && e.key == Makie.Keyboard.escape && e.action == Makie.Keyboard.press ||
            return Makie.Consume(false)
        close_context_menu!(v)
        return Makie.Consume(true)
    end
    Makie.on(ev.scroll; priority=MENU_PRIORITY) do _
        close_context_menu!(v)
        return Makie.Consume(false)
    end
end
