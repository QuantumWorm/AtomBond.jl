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

function build_panel!(pos, v::Viewer)
    vs = v.settings
    g = Makie.GridLayout(pos; width=280, valign=:top, tellheight=false, rowgap=4)
    row = 0
    nextrow() = (row += 1)
    heading(text) = Makie.Label(g[nextrow(), 1:3], text; font=:bold, halign=:left, fontsize=15,
                                padding=(0, 0, 0, row > 1 ? 6 : 0))
    button(pos, label) = with_click_feedback!(Makie.Button(pos; label, width=Makie.Relative(1),
                                                           fontsize=PANEL_FONTSIZE, padding=(6, 6, 5, 5)))
    function textbox(pos, placeholder, validator)
        t = Makie.Textbox(pos; placeholder, validator, width=Makie.Relative(1), fontsize=PANEL_FONTSIZE,
                          textpadding=(8, 8, 5, 5))
        sync_display!(t)
        push!(v.textboxes, t)
        return t
    end
    # Run `f(text)` with a text box's accepted value, or ask for one.
    function with_text(f, t, what)
        s, ok = box_text(t)
        isempty(s) ? (v.status[] = "Type $what first") :
        !ok ? (v.status[] = "Cannot use '$s': type $what") : f(s)
    end

    # --- Display style -------------------------------------------------------
    heading("Display style")
    menu = Makie.Menu(g[nextrow(), 1:3]; options=[(STYLE_LABELS[s], s) for s in DISPLAY_STYLES],
                      i_selected=findfirst(==(vs.style), DISPLAY_STYLES), fontsize=PANEL_FONTSIZE)
    Makie.on(menu.selection) do s
        vs.style = s
        redraw!(v)
    end

    # --- Show / hide -----------------------------------------------------------
    heading("Show")
    # Two columns of (switch, label) pairs.
    switches = Any[(label, x -> (setfield!(vs, field, x); redraw!(v)), getfield(vs, field))
                   for (label, field) in DISPLAY_TOGGLES]
    if v.triad !== nothing
        push!(switches, ("Axes", x -> (vs.show_triad = x; foreach(p -> p.visible = x, v.triad.plots)),
                         vs.show_triad))
    end
    tg = Makie.GridLayout(g[nextrow(), 1:3]; rowgap=6, colgap=6)
    for (k, (label, action, active)) in enumerate(switches)
        r, c = divrem(k - 1, 2) .+ (1, 0)
        t = Makie.Toggle(tg[r, 2c+1]; active)
        Makie.Label(tg[r, 2c+2], label; halign=:left, fontsize=PANEL_FONTSIZE)
        Makie.on(action, t.active)
    end

    # --- Sizes (labelled sliders, no heading) ---------------------------------
    sliders = Makie.SliderGrid(g[nextrow(), 1:3],
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

    # --- Camera ------------------------------------------------------------------
    heading("View along")
    vg = Makie.GridLayout(g[nextrow(), 1:3]; colgap=4)
    for (col, (label, which)) in enumerate(VIEW_BUTTONS)
        b = button(vg[1, col], label)
        Makie.on(_ -> view_along!(v, which), b.clicks)
    end
    r = nextrow()
    uvw = textbox(g[r, 1:2], "u v w  (e.g. 1 1 0)", is_int3)
    go = button(g[r, 3], "Go")
    Makie.on(go.clicks) do _
        with_text(uvw, "a direction u v w") do s
            dir = parse3(Int, s)
            if iszero(dir)
                v.status[] = "[0 0 0] is not a direction"
            else
                view_along!(v, dir)
                v.status[] = "Viewing along [$(join(dir, " "))]"
            end
        end
    end
    r = nextrow()
    reset = button(g[r, 1], "Reset")
    Makie.on(_ -> reset_camera!(v), reset.clicks)
    og = Makie.GridLayout(g[r, 2:3]; colgap=6)
    ortho = Makie.Toggle(og[1, 1]; active=vs.orthographic)
    Makie.Label(og[1, 2], "Orthographic"; halign=:left, fontsize=PANEL_FONTSIZE)
    Makie.on(ortho.active) do x
        vs.orthographic = x
        reframe!(v.lscene, v.display, vs)
    end

    # --- Cell ----------------------------------------------------------------------
    heading("Cell")
    r = nextrow()
    sc = textbox(g[r, 1], "n1 n2 n3", is_count3)
    apply = button(g[r, 2], "Supercell")
    surf = button(g[r, 3], "Surface…")
    Makie.on(_ -> surface_dialog(v), surf.clicks)
    Makie.on(apply.clicks) do _
        with_text(sc, "the repeats n1 n2 n3") do s
            n = parse3(Int, s)
            build!(v, "make $(join(n, "×")) supercell") do c
                any(iszero, n) && throw(ArgumentError("repeats must be ≥ 1"))
                supercell(c, n...)
            end
        end
    end
    r = nextrow()
    for (col, (label, f)) in enumerate([("Primitive", primitive_cell), ("Conventional", conventional_cell),
                                       ("Original", _ -> v.original)])
        b = button(g[r, col], label)
        Makie.on(_ -> build!(f, v, "go to $(lowercase(label)) cell"), b.clicks)
    end

    # --- Edit ----------------------------------------------------------------------
    heading("Edit selected atoms")
    r = nextrow()
    el = textbox(g[r, 1], "element", is_element)
    repl = button(g[r, 2], "Replace")
    Makie.on(_ -> with_text(s -> substitute_selected!(v, known_element(s)), el, "an element"), repl.clicks)
    del = button(g[r, 3], "Delete")
    Makie.on(_ -> delete_selected!(v), del.clicks)
    r = nextrow()
    shift = textbox(g[r, 1:2], "move by dx dy dz (Å)", is_real3)
    move = button(g[r, 3], "Move")
    Makie.on(_ -> with_text(s -> move_selected!(v, parse3(Float64, s)), shift, "dx dy dz"), move.clicks)
    r = nextrow()
    at = textbox(g[r, 1:2], "add at x y z (fractional)", is_real3)
    add = button(g[r, 3], "Add")
    Makie.on(add.clicks) do _
        with_text(el, "an element (box above)") do e
            with_text(at, "a position x y z") do s
                sym = known_element(e)
                edit!(c -> add_atom(c, sym, parse3(Float64, s)), v, "add $sym at ($(join(split(s), ", ")))")
            end
        end
    end
    r = nextrow()
    undo = button(g[r, 1], "Undo")
    Makie.on(_ -> undo!(v), undo.clicks)
    redo = button(g[r, 2], "Redo")
    Makie.on(_ -> redo!(v), redo.clicks)
    atoms_btn = button(g[r, 3], "Atoms…")   # edit element/position, replace a species everywhere
    Makie.on(_ -> atoms_dialog(v), atoms_btn.clicks)

    # --- Save ----------------------------------------------------------------------
    heading("File")
    r = nextrow()
    fg = Makie.GridLayout(g[r, 1:3]; colgap=6)
    # More in the File menu: recent files, Save as, Save image.
    newb = button(fg[1, 1], "New…")
    Makie.on(_ -> new_structure_dialog(v), newb.clicks)
    openb = button(fg[1, 2], "Open…")
    Makie.on(_ -> open_dialog(v), openb.clicks)
    saveb = button(fg[1, 3], "Save")
    Makie.on(_ -> save!(v), saveb.clicks)
    export_btn = button(fg[1, 4], "Export…")  # QE, CASTEP, VASP inputs, CIF, XYZ
    Makie.on(_ -> export_dialog(v), export_btn.clicks)

    Makie.colsize!(g, 1, Makie.Relative(1 / 3))
    return g
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
