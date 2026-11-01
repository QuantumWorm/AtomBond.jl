# Control panel on the right of the main window.

# (label, ViewSettings field) for the on/off switches that need a redraw.
# "Bond stubs" are the half-bonds to atoms outside the cell.
const DISPLAY_TOGGLES = [
    ("Bonds", :show_bonds),
    ("Bond stubs", :bond_stubs),
    ("Boundary atoms", :show_boundary_atoms),
    ("Unit cell", :show_cell),
    ("Atom labels", :show_labels),
]

const VIEW_BUTTONS = [("a", :a), ("b", :b), ("c", :c), ("a*", :astar), ("b*", :bstar), ("c*", :cstar)]

# Three numbers typed in a text box: separated by spaces, commas or
# semicolons; integers also by x / × ("2x2x2").
const REAL_SEP = r"[\s,;]+"
const INT_SEP = r"[\s,;x×*]+"

"""The three numbers in `s`, or `nothing`."""
function tryparse3(T::Type, s)
    w = split(strip(s), T <: Integer ? INT_SEP : REAL_SEP; keepempty=false)
    length(w) == 3 || return nothing
    x = tryparse.(T, w)
    return any(isnothing, x) ? nothing : T[x...]
end
function parse3(T, s)
    x = tryparse3(T, s)
    x === nothing && throw(ArgumentError("expected three numbers, got '$s'"))
    return x
end

# Text box validators.
is_int3(s) = tryparse3(Int, s) !== nothing
is_count3(s) = (n = tryparse3(Int, s)) !== nothing && all(>=(0), n)
is_real3(s) = tryparse3(Float64, s) !== nothing
is_element(s) = try
    known_element(s); true
catch
    false
end

"""
Show a value stored in a text box from code (not by typing + Enter), so
[`box_text`](@ref), which reads what is shown, sees it too.
"""
function sync_display!(t)
    Makie.on(t.stored_string) do s
        s === nothing || t.focused[] || t.displayed_string[] == s || (t.displayed_string[] = s)
    end
    return t
end

"""
    box_text(textbox) -> (text, ok)

What a text box shows, so typed text counts even before Enter is pressed,
and whether the box's validator accepts it. An empty box gives ("", false).
"""
function box_text(t)
    s = strip(something(t.displayed_string[], ""))
    isempty(s) && return s, false
    val = t.validator[]
    ok = val isa Function ? val(s) : val isa Regex ? occursin(val, s) : true
    return s, ok
end

# Compact controls so the whole panel fits a ~960 px tall window.
const PANEL_FONTSIZE = 14

"""Make a structure-changing step from the Cell section; these re-frame the camera."""
build!(f, v::Viewer, what::AbstractString) = edit!(f, v, what; reframe=true)

"""The view controls being built: their three-column layout, the last row used and the viewer."""
mutable struct Panel
    g::Makie.GridLayout
    row::Int
    v::Viewer
end

nextrow!(p::Panel) = (p.row += 1)

heading!(p::Panel, text) = Makie.Label(p.g[nextrow!(p), 1:3], text; font=:bold, halign=:left, fontsize=15,
                                       padding=(0, 0, 0, p.row > 1 ? 6 : 0))

panel_button(pos, label) = with_click_feedback!(Makie.Button(pos; label, width=Makie.Relative(1),
                                                             fontsize=PANEL_FONTSIZE, padding=(6, 6, 5, 5)))

function panel_textbox!(p::Panel, pos, placeholder, validator)
    t = Makie.Textbox(pos; placeholder, validator, width=Makie.Relative(1), fontsize=PANEL_FONTSIZE,
                      textpadding=(8, 8, 5, 5))
    sync_display!(t)
    push!(p.v.textboxes, t)
    return t
end

"""Run `f(text)` with a text box's accepted value, or ask for one in the status line."""
function with_text(f, v::Viewer, t, what)
    s, ok = box_text(t)
    isempty(s) ? (v.status[] = "Type $what first") :
    !ok ? (v.status[] = "Cannot use '$s': type $what") : f(s)
end

"""The view controls: display style, what is shown, sizes, view direction, cell, editing and file."""
function build_panel!(pos, v::Viewer)
    g = Makie.GridLayout(pos; width=280, valign=:top, tellheight=false, rowgap=4)
    p = Panel(g, 0, v)
    style_section!(p)
    show_section!(p)
    size_section!(p)
    view_section!(p)
    cell_section!(p)
    edit_section!(p)
    file_section!(p)
    Makie.colsize!(g, 1, Makie.Relative(1 / 3))
    return g
end

function style_section!(p::Panel)
    v, vs = p.v, p.v.settings
    heading!(p, "Display style")
    menu = Makie.Menu(p.g[nextrow!(p), 1:3]; options=[(STYLE_LABELS[s], s) for s in DISPLAY_STYLES],
                      i_selected=findfirst(==(vs.style), DISPLAY_STYLES), fontsize=PANEL_FONTSIZE)
    Makie.on(menu.selection) do s
        vs.style = s
        redraw!(v)
    end
end

"""On/off switches in two columns of (switch, label) pairs."""
function show_section!(p::Panel)
    v, vs = p.v, p.v.settings
    heading!(p, "Show")
    switches = Any[(label, x -> (setfield!(vs, field, x); redraw!(v)), getfield(vs, field))
                   for (label, field) in DISPLAY_TOGGLES]
    if v.triad !== nothing
        push!(switches, ("Axes", x -> (vs.show_triad = x; foreach(t -> t.visible = x, v.triad.plots)),
                         vs.show_triad))
    end
    tg = Makie.GridLayout(p.g[nextrow!(p), 1:3]; rowgap=6, colgap=6)
    for (k, (label, action, active)) in enumerate(switches)
        r, c = divrem(k - 1, 2) .+ (1, 0)
        t = Makie.Toggle(tg[r, 2c+1]; active)
        Makie.Label(tg[r, 2c+2], label; halign=:left, fontsize=PANEL_FONTSIZE)
        Makie.on(action, t.active)
    end
end

"""Atom and bond size sliders (no heading)."""
function size_section!(p::Panel)
    v, vs = p.v, p.v.settings
    sliders = Makie.SliderGrid(p.g[nextrow!(p), 1:3],
        (label="Atoms", range=0.2:0.05:1.0, startvalue=vs.ball_scale, format="{:.2f}"),
        (label="Bonds", range=0.04:0.01:0.30, startvalue=vs.bond_radius, format="{:.2f} Å");
        tellheight=true)
    Makie.on(sliders.sliders[1].value) do x
        vs.ball_scale = x
        vs.cpk_scale = x / 0.4      # same slider scales CPK spheres around 1.0
        redraw!(v)
    end
    Makie.on(sliders.sliders[2].value) do x
        vs.bond_radius = x
        redraw!(v)
    end
end

"""View along a lattice or reciprocal axis or [uvw], turn the view in the screen plane, reset, orthographic."""
function view_section!(p::Panel)
    v, vs, g = p.v, p.v.settings, p.g
    heading!(p, "View along")
    vg = Makie.GridLayout(g[nextrow!(p), 1:3]; colgap=4)
    for (col, (label, which)) in enumerate(VIEW_BUTTONS)
        b = panel_button(vg[1, col], label)
        Makie.on(_ -> view_along!(v, which), b.clicks)
    end
    r = nextrow!(p)
    uvw = panel_textbox!(p, g[r, 1:2], "u v w  (e.g. 1 1 0)", is_int3)
    go = panel_button(g[r, 3], "Go")
    Makie.on(go.clicks) do _
        with_text(v, uvw, "a direction u v w") do s
            dir = parse3(Int, s)
            if iszero(dir)
                v.status[] = "[0 0 0] is not a direction"
            else
                view_along!(v, dir)
                v.status[] = "Viewing along [$(join(dir, " "))]"
            end
        end
    end
    # Turn the picture in the screen plane (the line of sight stays).
    r = nextrow!(p)
    ccw = panel_button(g[r, 1], "↺ 90°")
    cw = panel_button(g[r, 2], "↻ 90°")
    Makie.on(_ -> (rotate_view!(v, 90); v.status[] = "View turned 90° counter-clockwise"), ccw.clicks)
    Makie.on(_ -> (rotate_view!(v, -90); v.status[] = "View turned 90° clockwise"), cw.clicks)
    rg = Makie.GridLayout(g[r, 3]; colgap=4)
    angle = panel_textbox!(p, rg[1, 1], "°", s -> isnum(s))
    turn = panel_button(rg[1, 2], "Turn")
    Makie.on(turn.clicks) do _
        with_text(v, angle, "an angle in degrees (+ counter-clockwise)") do s
            x = parse(Float64, strip(s))
            rotate_view!(v, x)
            v.status[] = "View turned $(x)° " * (x >= 0 ? "counter-clockwise" : "clockwise")
        end
    end
    r = nextrow!(p)
    reset = panel_button(g[r, 1], "Reset")
    Makie.on(_ -> reset_camera!(v), reset.clicks)
    og = Makie.GridLayout(g[r, 2:3]; colgap=6)
    ortho = Makie.Toggle(og[1, 1]; active=vs.orthographic)
    Makie.Label(og[1, 2], "Orthographic"; halign=:left, fontsize=PANEL_FONTSIZE)
    Makie.on(ortho.active) do x
        vs.orthographic = x
        reframe!(v.lscene, v.display, vs)
    end
end

"""Supercell, surface, primitive / conventional / original cell."""
function cell_section!(p::Panel)
    v, g = p.v, p.g
    heading!(p, "Cell")
    r = nextrow!(p)
    sc = panel_textbox!(p, g[r, 1], "n1 n2 n3", is_count3)
    apply = panel_button(g[r, 2], "Supercell")
    surf = panel_button(g[r, 3], "Surface…")
    Makie.on(_ -> surface_dialog(v), surf.clicks)
    Makie.on(apply.clicks) do _
        with_text(v, sc, "the repeats n1 n2 n3") do s
            n = parse3(Int, s)
            build!(v, "make $(join(n, "×")) supercell") do c
                any(iszero, n) && throw(ArgumentError("repeats must be ≥ 1"))
                supercell(c, n...)
            end
        end
    end
    r = nextrow!(p)
    for (col, (label, f)) in enumerate([("Primitive", primitive_cell), ("Conventional", conventional_cell),
                                       ("Original", _ -> v.original)])
        b = panel_button(g[r, col], label)
        Makie.on(_ -> build!(f, v, "go to $(lowercase(label)) cell"), b.clicks)
    end
end

"""Replace, delete, move or add atoms; undo and redo; the Edit atoms tool."""
function edit_section!(p::Panel)
    v, g = p.v, p.g
    heading!(p, "Edit selected atoms")
    r = nextrow!(p)
    el = panel_textbox!(p, g[r, 1], "element", is_element)
    repl = panel_button(g[r, 2], "Replace")
    Makie.on(_ -> with_text(s -> substitute_selected!(v, known_element(s)), v, el, "an element"), repl.clicks)
    del = panel_button(g[r, 3], "Delete")
    Makie.on(_ -> delete_selected!(v), del.clicks)
    r = nextrow!(p)
    shift = panel_textbox!(p, g[r, 1:2], "move by dx dy dz (Å)", is_real3)
    move = panel_button(g[r, 3], "Move")
    Makie.on(_ -> with_text(s -> move_selected!(v, parse3(Float64, s)), v, shift, "dx dy dz"), move.clicks)
    r = nextrow!(p)
    at = panel_textbox!(p, g[r, 1:2], "add at x y z (fractional)", is_real3)
    add = panel_button(g[r, 3], "Add")
    Makie.on(add.clicks) do _
        with_text(v, el, "an element (box above)") do e
            with_text(v, at, "a position x y z") do s
                sym = known_element(e)
                edit!(c -> add_atom(c, sym, parse3(Float64, s)), v, "add $sym at ($(join(split(s), ", ")))")
            end
        end
    end
    r = nextrow!(p)
    undo = panel_button(g[r, 1], "Undo")
    Makie.on(_ -> undo!(v), undo.clicks)
    redo = panel_button(g[r, 2], "Redo")
    Makie.on(_ -> redo!(v), redo.clicks)
    atoms_btn = panel_button(g[r, 3], "Atoms…")   # edit element/position, replace a species everywhere
    Makie.on(_ -> atoms_dialog(v), atoms_btn.clicks)
end

"""New, Open, Save, Export (more in the File menu: recent files, Save as, Save image)."""
function file_section!(p::Panel)
    v = p.v
    heading!(p, "File")
    fg = Makie.GridLayout(p.g[nextrow!(p), 1:3]; colgap=6)
    for (col, (label, action)) in enumerate([("New…", new_structure_dialog), ("Open…", open_dialog),
                                             ("Save", save!), ("Export…", export_dialog)])
        b = panel_button(fg[1, col], label)
        Makie.on(_ -> action(v), b.clicks)
    end
end

# --- Docked tools -------------------------------------------------------------------
#
# The right-hand column (figure rows 1–3, column 2) holds either the view
# controls above or one tool (Supercell, Shift atoms, ...). Makie blocks
# cannot be hidden from the mouse, so the column's contents are deleted and
# rebuilt when it switches.

"""Delete every block of a layout (also those inside blocks like SliderGrid) and the layout itself."""
function delete_layout!(g::Makie.GridLayout)
    for c in copy(g.content)
        x = c.content
        if x isa Makie.GridLayout
            delete_layout!(x)
        else
            hasproperty(x, :layout) && x.layout isa Makie.GridLayout && delete_layout!(x.layout)
            delete!(x)
        end
    end
    gc = Makie.GridLayoutBase.gridcontent(g)
    gc === nothing || Makie.GridLayoutBase.remove_from_gridlayout!(gc)
    return nothing
end

"""Forget text boxes that are being deleted, so shortcuts no longer check them."""
forget_textboxes!(v::Viewer, boxes) = filter!(t -> !any(b -> b === t, boxes), v.textboxes)

"""Show the view controls in the right-hand column (if the viewer has them)."""
function show_panel!(v::Viewer)
    haskey(v.dialogs, :panel) && return nothing
    n = length(v.textboxes)
    g = build_panel!(v.figure[1:3, 2], v)
    v.dialogs[:panel] = (; layout=g, boxes=v.textboxes[n+1:end])
    return g
end

function hide_panel!(v::Viewer)
    p = pop!(v.dialogs, :panel, nothing)
    p === nothing && return nothing
    forget_textboxes!(v, p.boxes)
    delete_layout!(p.layout)
end

"""
Clear the right-hand column and start a tool there: a title with a Close
button (and the tab it works on, `subject`), then an empty Form `width` px
wide for the tool's controls.
"""
function dock_tool!(v::Viewer, title::AbstractString, width::Real; subject=nothing)
    close_tool!(v; restore=false)
    hide_panel!(v)
    g = Makie.GridLayout(v.figure[1:3, 2]; width, valign=:top, tellheight=false, rowgap=6)
    top = Makie.GridLayout(g[1, 1:2]; colgap=8)
    Makie.Label(top[1, 1], title; font=:bold, fontsize=17, halign=:left, tellwidth=false)
    close_btn = with_click_feedback!(Makie.Button(top[1, 2]; label="Close", fontsize=PANEL_FONTSIZE,
                                                  padding=(10, 10, 4, 4)))
    Makie.on(_ -> close_tool!(v), close_btn.clicks)
    subject === nothing && return Form(g, 1)
    Makie.Label(g[2, 1:2], "on $subject"; halign=:left, fontsize=13, color=:gray40)
    return Form(g, 2)
end

"""
    close_tool!(viewer; restore=true) -> Bool

Close the tool open on the right; `restore` brings the view controls back.
Returns whether a tool was open.
"""
function close_tool!(v::Viewer; restore::Bool=true)
    t = pop!(v.dialogs, :tool, nothing)
    if t !== nothing
        forget_textboxes!(v, t.boxes)
        delete_layout!(t.layout)
    end
    restore && get(v.dialogs, :has_panel, false) && show_panel!(v)
    return t !== nothing
end

"""Key of the tool open on the right (:supercell, :shift, ...), or `nothing`."""
open_tool(v::Viewer) = (t = get(v.dialogs, :tool, nothing)) === nothing ? nothing : t.key
