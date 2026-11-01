# Tools of the main viewer. Each opens in the right-hand column of the main
# window, in place of the view controls, with a Close button that brings the
# view controls back (one tool at a time; see `dock_tool!` in panel.jl):
#
# - Build surface: cleave slabs from, and build superlattices on, the current
#   structure. Results are undoable edits (Ctrl+Z in the main window).
# - Oxygen bridge: put an O (or another atom) into a selected bond.
# - Export: write input files for VASP / Quantum ESPRESSO / CASTEP, or CIF / XYZ.
# - Open structure: a file browser (filebrowser.jl).
# - Edit atoms: element / position of selected atoms, replace a species everywhere.
# - Supercell: repeats n1 n2 n3 or new cell vectors as integer combinations.
# - Shift atoms: translate all or the selected atoms; wrap them into the cell.
# - Rotate atoms: turn all or the selected atoms about an axis and a centre.
# - Structure library: open a known crystal or molecule (io/library.jl).
# - Cell and constraints: cell data, periodic box for molecules, lattice and
#   atom constraints for geometry optimization.
#
# Messages appear in the main window's status line.

const DIALOG_WIDTH = 440

"""Two-column form (label | control) used by the dialogs; `boxes` are its text boxes."""
mutable struct Form
    g::Makie.GridLayout
    row::Int
    boxes::Vector{Any}
end
Form(g, row) = Form(g, row, Any[])

nextrow!(f::Form) = (f.row += 1)

heading!(f::Form, text) = Makie.Label(f.g[nextrow!(f), 1:2], text; font=:bold, halign=:left, fontsize=15,
                                      padding=(0, 0, 0, f.row > 1 ? 10 : 0))

note!(f::Form, text) = Makie.Label(f.g[nextrow!(f), 1:2], text; halign=:left, fontsize=12,
                                   color=:gray40, justification=:left)

function field!(f::Form, label, default, validator)
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], label; halign=:left, fontsize=PANEL_FONTSIZE)
    t = sync_display!(Makie.Textbox(f.g[r, 2]; stored_string=default, validator, width=Makie.Relative(1),
                                    fontsize=PANEL_FONTSIZE, textpadding=(8, 8, 5, 5)))
    push!(f.boxes, t)
    return t
end

button!(f::Form, label) = with_click_feedback!(Makie.Button(f.g[nextrow!(f), 1:2]; label, width=Makie.Relative(1),
                                                            fontsize=PANEL_FONTSIZE, padding=(6, 6, 5, 5)))

isnum(s) = tryparse(Float64, strip(s)) !== nothing

"""A text box's accepted value, or `nothing` after asking for it in the status line."""
function form_value(v::Viewer, t, what)
    s, ok = box_text(t)
    isempty(s) ? (v.status[] = "Enter $what"; nothing) :
    !ok ? (v.status[] = "Cannot use '$s' for $what"; nothing) : s
end

"""
Start a tool: an empty Form docked in the main window's right-hand column
(replacing the view controls or another tool), or `nothing` (with a
message) if this tool is already open.
"""
function new_dialog(v::Viewer, key::Symbol, title::AbstractString; width::Integer=DIALOG_WIDTH)
    tool = get(v.dialogs, :tool, nothing)
    if tool !== nothing && tool.key === key
        v.status[] = "$title is already open on the right"
        return nothing
    end
    subject = haskey(STRUCTURE_TOOLS, key) ? tab_label(v, v.active) : nothing   # the tab it works on
    return dock_tool!(v, title, width - 32; subject)
end

"""
Break `s` into lines of at most `n` characters at spaces (existing line
breaks stay; a word longer than `n` gets a line of its own).
"""
function wrap_text(s::AbstractString, n::Integer)
    lines = String[]
    for para in split(s, '\n')
        line = ""
        for word in split(para, ' ')
            if isempty(line) || length(line) + 1 + length(word) <= n
                line = isempty(line) ? String(word) : line * " " * word
            else
                push!(lines, line)
                line = String(word)
            end
        end
        push!(lines, line)
    end
    return join(lines, '\n')
end

"""
Wrap the text of a form's labels at the width of their columns (the form is
`width` px wide, column 1 takes `labelwidth` of it), now and whenever the
text changes, so long lines do not widen the column. Makie's own word_wrap,
or a fixed label width, hangs the layout here.
"""
function wrap_labels!(f::Form, width::Real, labelwidth::Real)
    for c in f.g.content
        x = c.content
        x isa Makie.Label || continue
        cols = c.span.cols
        w = cols == 1:2 ? width : first(cols) == 1 ? labelwidth * width : (1 - labelwidth) * width - 10
        n = floor(Int, w / (0.55 * x.fontsize[]))          # characters per line, a little conservatively
        rewrap(s) = s isa AbstractString && (t = wrap_text(s, n)) != s && (x.text[] = t)
        x.justification[] = :left
        Makie.on(rewrap, x.text)
        rewrap(x.text[])
    end
end

"""Finish a tool's form; it becomes the open tool of the viewer."""
function finish_dialog!(v::Viewer, key::Symbol, title::AbstractString, f::Form; labelwidth::Real=0.4)
    Makie.colsize!(f.g, 1, Makie.Relative(labelwidth))
    Makie.colsize!(f.g, 2, Makie.Auto(false))     # take the remaining width, not the widest item
    wrap_labels!(f, f.g.width[], labelwidth)      # long lines wrap instead of widening the column
    append!(v.textboxes, f.boxes)                 # keyboard shortcuts pause while typing in them
    v.dialogs[:tool] = (; key, title, layout=f.g, boxes=f.boxes)
end

# ---------------------------------------------------------------------------

"""
    surface_dialog(viewer) -> NamedTuple of the tool's controls

Open the Build surface tool for `viewer`.
"""
function surface_dialog(v::Viewer)
    f = new_dialog(v, :surface, "Build surface")
    f === nothing && return nothing

    heading!(f, "Cleave surface (current structure)")
    hkl = field!(f, "Plane  h k l", "0 0 1", s -> is_int3(s) && !iszero(parse3(Int, s)))
    layers = field!(f, "Layers", "4", s -> something(tryparse(Int, strip(s)), 0) >= 1)
    vacuum = field!(f, "Vacuum (Å)", "15", s -> isnum(s) && parse(Float64, s) >= 0)
    shift = field!(f, "Termination shift (0–1)", "auto", s -> isnum(s) || strip(lowercase(s)) == "auto")
    note!(f, "One layer = one repeat of the surface cell along the normal.\n" *
             "Shift \"auto\" cuts where the fewest bonds break; a number\n" *
             "(fraction of a layer) moves the cut to another atomic plane.")
    build_slab = button!(f, "Build slab")

    heading!(f, "Superlattice  A / B")
    note!(f, "A = current structure, B = a file; stacked along the plane above.")
    bfile = field!(f, "B file", "", s -> isfile(expanduser(strip(s))))
    ab = field!(f, "Layers  A B", "4 4", s -> occursin(r"^\s*\d+\s+\d+\s*$", s) && all(>(0), parse.(Int, split(s))))
    r = nextrow!(f)
    strain = Makie.Toggle(f.g[r, 1]; active=true, halign=:left)
    Makie.Label(f.g[r, 2], "Strain B in plane to A"; halign=:left, fontsize=PANEL_FONTSIZE)
    build_sl = button!(f, "Build superlattice")


    Makie.on(build_slab.clicks) do _
        vals = (form_value(v, hkl, "the plane h k l"), form_value(v, layers, "the number of layers"),
                form_value(v, vacuum, "the vacuum"), form_value(v, shift, "the termination shift"))
        any(isnothing, vals) && return
        p = parse3(Int, vals[1])
        n, vac = parse(Int, vals[2]), parse(Float64, vals[3])
        sh = lowercase(vals[4]) == "auto" ? :auto : parse(Float64, vals[4])
        edit!(c -> slab(c, p; layers=n, vacuum=vac, shift=sh), v,
              "cleave ($(join(p, " "))) slab, $n layers, $(vac) Å vacuum"; reframe=true)
    end

    Makie.on(build_sl.clicks) do _
        vals = (form_value(v, hkl, "the plane h k l"), form_value(v, bfile, "a file for B"),
                form_value(v, ab, "layers A B"))
        any(isnothing, vals) && return
        p = parse3(Int, vals[1])
        nA, nB = parse.(Int, split(vals[3]))
        path = expanduser(vals[2])
        edit!(v, "build ($(join(p, " "))) superlattice with $(basename(path))"; reframe=true) do c
            superlattice(c, load_structure(path), p; layers=(nA, nB), strain=strain.active[])
        end
    end

    finish_dialog!(v, :surface, "Build surface", f)
    return (; hkl, layers, vacuum, shift, build_slab, bfile, ab, strain, build_sl)
end

# ---------------------------------------------------------------------------

"""
    export_dialog(viewer) -> NamedTuple of the tool's controls

Open the Export tool: choose a format, folder and name, k-point spacing
(the resulting mesh for the current structure is shown) and, for Quantum
ESPRESSO, ABINIT and DFTK, the cutoff and pseudopotentials.
"""
function export_dialog(v::Viewer)
    f = new_dialog(v, :export, "Export")
    f === nothing && return nothing

    heading!(f, "Export current structure")
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], "Format"; halign=:left, fontsize=PANEL_FONTSIZE)
    format = Makie.Menu(f.g[r, 2]; options=EXPORT_FORMATS, i_selected=1, fontsize=PANEL_FONTSIZE,
                        width=Makie.Relative(1))
    folder = field!(f, "Folder", pwd(), s -> !isempty(strip(s)))
    name = field!(f, "Name", formula(v.crystal), s -> occursin(r"^[\w.+-]+$", strip(s)))
    note!(f, "Files are Folder/Name.<ext>; VASP writes a folder Folder/Name/.")

    heading!(f, "k-points")
    kspacing = field!(f, "Spacing (1/Å, incl. 2π)", "0.25", s -> isnum(s) && parse(Float64, s) > 0)
    # Mesh for the current structure, updated as the spacing or structure changes.
    mesh = Makie.Observable("")
    function update_mesh(_...)
        s = kspacing.stored_string[]
        mesh[] = s === nothing || !isnum(s) ? "" :
                 "Mesh for this structure: " * join(kpoint_grid(v.crystal; kspacing=parse(Float64, s)), " × ")
    end
    Makie.on(update_mesh, kspacing.stored_string)
    Makie.on(update_mesh, v.header)          # the header changes whenever the structure does
    update_mesh()
    Makie.Label(f.g[nextrow!(f), 1:2], mesh; halign=:left, fontsize=PANEL_FONTSIZE)

    heading!(f, "Plane waves (Quantum ESPRESSO, ABINIT, DFTK)")
    ecut = field!(f, "Cutoff (Ry)", "40", s -> isnum(s) && parse(Float64, s) > 0)
    pseudo = field!(f, "Pseudopotentials", "{El}.UPF", s -> occursin("{El}", s) || occursin(r"^[\w.+-]+$", s))
    note!(f, "{El} is replaced by each element symbol; DFTK also takes\na PseudoPotentialData family name. ABINIT and DFTK\nget the cutoff in Ha (half the Ry value).")
    # The usual pseudopotential files follow the code.
    defaults = ("{El}.UPF", "{El}.psp8", DFTK_PSEUDO_FAMILY)
    Makie.on(format.selection) do fmt
        t = shown_text(pseudo)
        t in defaults || return
        fmt === :abinit && settext!(pseudo, "{El}.psp8")
        fmt === :qe && settext!(pseudo, "{El}.UPF")
        fmt === :dftk && settext!(pseudo, DFTK_PSEUDO_FAMILY)
    end

    go = button!(f, "Export")
    Makie.on(go.clicks) do _
        vals = (form_value(v, folder, "a folder"), form_value(v, name, "a name"),
                form_value(v, kspacing, "a k-point spacing"), form_value(v, ecut, "ecutwfc"),
                form_value(v, pseudo, "a pseudopotential pattern"))
        any(isnothing, vals) && return
        try
            paths = export_inputs(v.crystal, format.selection[], expanduser(vals[1]), vals[2];
                                  kspacing=parse(Float64, vals[3]), ecutwfc=parse(Float64, vals[4]),
                                  pseudo=vals[5])
            v.status[] = "Wrote " * join(basename.(paths), ", ") * " in " * dirname(first(paths))
        catch err
            err isa InterruptException && rethrow()
            v.status[] = "Could not export: " * sprint(showerror, err)
        end
    end

    finish_dialog!(v, :export, "Export", f)
    return (; format, folder, name, kspacing, mesh, ecut, pseudo, go)
end

# ---------------------------------------------------------------------------

const COORD_MODES = ["Fractional", "Cartesian (Å)"]

"""Coordinates of the shown copy (index, image) of an atom: fractional or Cartesian Å."""
function atom_coordinates(c::Crystal, index::Integer, image, cartesian::Bool)
    f = c.frac[index] + Vec3(image)
    return cartesian ? frac_to_cart(c.lattice, f) : f
end

fmt3(x) = join((@sprintf("%.6f", t) for t in x), "  ")

# Text box contents as shown (typed text counts even before Enter is pressed).
shown_text(t) = strip(something(t.displayed_string[], ""))

function settext!(t, s::AbstractString)
    t.displayed_string[] = s
    t.stored_string[] = s
end

"""
    atoms_dialog(viewer) -> NamedTuple of the tool's controls

Edit atoms tool. Selected atoms: change the element (all selected atoms)
and, for a single atom, set its position (fractional or Cartesian Å; for an
image atom, ′, the position of that copy). Replace species: change every atom
of one element into another. All changes are undoable.
"""
function atoms_dialog(v::Viewer)
    f = new_dialog(v, :atoms, "Edit atoms")
    f === nothing && return nothing

    heading!(f, "Selected atoms")
    who = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], who; halign=:left, fontsize=PANEL_FONTSIZE)
    element = field!(f, "Element", "", s -> isempty(strip(s)) || is_element(s))
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], "Coordinates"; halign=:left, fontsize=PANEL_FONTSIZE)
    mode = Makie.Menu(f.g[r, 2]; options=COORD_MODES, i_selected=1, fontsize=PANEL_FONTSIZE, width=Makie.Relative(1))
    coords = field!(f, "x  y  z", "", s -> isempty(strip(s)) || is_real3(s))
    note!(f, "One atom: change its element and/or position. Several\n" *
             "atoms: the element is applied to all of them.")
    apply = button!(f, "Apply")

    heading!(f, "Replace a species everywhere")
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], "Replace"; halign=:left, fontsize=PANEL_FONTSIZE)
    from = Makie.Menu(f.g[r, 2]; options=string.(species_list(v.crystal)), i_selected=1,
                      fontsize=PANEL_FONTSIZE, width=Makie.Relative(1))
    to = field!(f, "With", "", s -> isempty(strip(s)) || is_element(s))
    count_text = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], count_text; halign=:left, fontsize=12, color=:gray40)
    replace_all = button!(f, "Replace all")

    function show_selection(_...)
        c = v.crystal
        idx = selected_indices(v)
        if isempty(idx)
            who[] = "None — click atoms in the main window"
            settext!(element, ""); settext!(coords, "")
        elseif length(idx) == 1
            i, img = first(v.selection)
            who[] = atom_label(c, i) * (iszero(img) ? "" : "′  (copy in a neighbouring cell)")
            settext!(element, string(c.species[i]))
            settext!(coords, fmt3(atom_coordinates(c, i, img, mode.i_selected[] == 2)))
        else
            spc = c.species[idx]
            who[] = "$(length(idx)) atoms: " * join(("$s×$(count(==(s), spc))" for s in unique(spc)), ", ")
            settext!(element, allequal(spc) ? string(spc[1]) : "")
            settext!(coords, "")
        end
        # Offer the species of the first selected atom for replacement.
        isempty(idx) || (k = findfirst(==(string(c.species[idx[1]])), from.options[])) === nothing ||
            (from.i_selected[] = k)
    end
    function show_species(_...)
        opts = string.(species_list(v.crystal))
        if opts != from.options[]
            from.options[] = opts
            from.i_selected[] = 1
        end
    end
    function show_count(_...)
        s = from.selection[]
        count_text[] = s === nothing ? "" : plural(count(==(Symbol(s)), v.crystal.species), "atom") * " of $s"
    end
    Makie.on(_ -> (show_species(); show_count()), v.header)       # structure changed
    Makie.on(show_selection, v.info)                              # selection changed
    Makie.on(show_selection, mode.i_selected)
    Makie.on(show_count, from.selection)
    show_species(); show_selection(); show_count()

    Makie.on(apply.clicks) do _
        idx = selected_indices(v)
        isempty(idx) && (v.status[] = "Select atoms first (click / Shift+click)"; return)
        es = shown_text(element)
        el = isempty(es) ? nothing : is_element(es) ? known_element(es) : (v.status[] = "Unknown element '$es'"; return)
        if length(idx) > 1
            el === nothing && (v.status[] = "Enter an element for the selected atoms"; return)
            return substitute_selected!(v, el)
        end
        c = v.crystal
        i, img = first(v.selection)
        cart = mode.i_selected[] == 2
        cs = shown_text(coords)
        is_real3(cs) || (v.status[] = "Enter the position as x y z"; return)
        p = Vec3(parse3(Float64, cs))
        moved = norm(p - atom_coordinates(c, i, img, cart)) > 1e-6
        changed = el !== nothing && el != c.species[i]
        moved || changed || (v.status[] = "Nothing to change"; return)
        # Coordinates are those of the copy shown; store the home-cell position.
        fhome = (cart ? cart_to_frac(c.lattice, p) : p) - Vec3(img)
        what = "edit " * atom_label(c, i) * (changed ? " → $el" : "") * (moved ? " at ($cs)" : "")
        edit!(v, what; keep_selection=true) do c
            c = changed ? substitute(c, [i], el) : c
            moved ? set_position(c, i, fhome) : c
        end
    end

    Makie.on(replace_all.clicks) do _
        src = from.selection[]
        ts = shown_text(to)
        (src === nothing || isempty(ts)) && (v.status[] = "Choose a species and the element to put in its place"; return)
        is_element(ts) || (v.status[] = "Unknown element '$ts'"; return)
        dst = known_element(ts)
        n = count(==(Symbol(src)), v.crystal.species)
        edit!(c -> replace_species(c, src, dst), v, "replace all $src with $dst ($(plural(n, "atom")))";
              keep_selection=true)
    end

    finish_dialog!(v, :atoms, "Edit atoms", f)
    return (; who, element, mode, coords, apply, from, to, count_text, replace_all)
end

# ---------------------------------------------------------------------------

# (name, a′, b′, c′) in units of a, b, c; c′ of the in-plane presets is set
# from the third repeat.
const SUPERCELL_PRESETS = [
    ("Repeat n1 × n2 × n3", nothing),
    ("√2×√2 R45° in ab (square)", ((1, 1, 0), (-1, 1, 0))),
    ("√3×√3 R30° in ab (hexagonal)", ((2, 1, 0), (-1, 1, 0))),
    ("Rectangular 1×√3 in ab (hexagonal)", ((1, 0, 0), (1, 2, 0))),
]

"""The supercell matrix P (columns a′, b′, c′) from three "x y z" texts, or `nothing`."""
function supercell_matrix(texts)
    cols = [tryparse3(Int, t) for t in texts]
    any(isnothing, cols) && return nothing
    return hcat(cols...)
end

"""One line about the supercell `P` of `c`: atoms and cell, or what is wrong."""
function supercell_preview(c::Crystal, P)
    P === nothing && return "Enter a′ b′ c′ as three integers each"
    n = round(Int, det(P))
    n > 0 || return "The new vectors must be right-handed and independent (det = $n)"
    a, b, cc, α, β, γ = lattice_parameters(Lattice(c.lattice.matrix * P))
    return @sprintf("%d× the cell: %d atoms\na b c = %.3f %.3f %.3f Å   α β γ = %.1f %.1f %.1f°",
                    n, n * natoms(c), a, b, cc, α, β, γ)
end

"""Two non-negative integers n m, not both zero (nanotube indices), or `nothing`."""
function nanotube_indices(s::AbstractString)
    t = split(replace(s, ',' => ' '))
    length(t) == 2 || return nothing
    nm = tryparse.(Int, t)
    any(isnothing, nm) && return nothing
    n, m = nm
    return n >= 0 && m >= 0 && n + m > 0 ? (n, m) : nothing
end

"""Preview text of the nanotube tool: kind, diameter, period and number of atoms, or what to fix."""
function nanotube_preview(nm, bond, periods, vacuum, els)
    nm === nothing && return "Enter the indices as two integers n m (n, m ≥ 0), e.g. 6 5"
    bond === nothing && return "Enter the bond length in Å"
    periods === nothing && return "Enter the number of periods (1 or more)"
    vacuum === nothing && return "Enter the vacuum in Å"
    els === nothing && return "Enter one element (C) or two, one per sublattice (B N)"
    n, m = nm
    _, _, R, T = nanotube_geometry(n, m, bond)
    shape, metallic = nanotube_kind(n, m)
    natom = periods * 4 * (n^2 + n * m + m^2) ÷ gcd(2m + n, 2n + m)
    elec = els[1] == els[2] == :C ? (metallic ? ", metallic" : ", semiconducting") : ""
    return @sprintf("(%d,%d) %s%s\ndiameter %.3f Å, period %.3f Å along c\n%d atoms in the cell (%d period%s)",
                    n, m, shape, elec, 2R, T, natom, periods, periods == 1 ? "" : "s")
end

"""
    nanotube_dialog(viewer) -> NamedTuple of the tool's controls

Nanotube tool: a single-wall (n, m) nanotube (carbon by default, or two
elements such as B N), periodic along c, opened in a new tab.
"""
function nanotube_dialog(v::Viewer)
    f = new_dialog(v, :nanotube, "Nanotube")
    f === nothing && return nothing
    heading!(f, "Single-wall (n, m) nanotube")
    indices = field!(f, "Indices n m", "6 5", s -> nanotube_indices(s) !== nothing)
    periods = field!(f, "Periods along c", "1", s -> (k = tryparse(Int, strip(s)); k !== nothing && k >= 1))
    bond = field!(f, "Bond length (Å)", @sprintf("%.4f", GRAPHITE_CC), s -> isnum(s) && parse(Float64, s) > 0)
    vacuum = field!(f, "Vacuum (Å)", "10", s -> isnum(s) && parse(Float64, s) > 0)
    elements = field!(f, "Elements", "C", s -> tube_elements(s) !== nothing)
    preview = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], preview; halign=:left, fontsize=13, color=:gray25, justification=:left)
    note!(f, "n = m armchair, m = 0 zigzag, else chiral.\n" *
             "Carbon: metallic when n − m is a multiple of 3.\n" *
             "Bond length: graphite's C–C by default.\n" *
             "An ideal rolled sheet: relax it before use.")
    create = button!(f, "Create nanotube")
    num(t) = (x = tryparse(Float64, strip(shown_text(t))); x === nothing || x <= 0 ? nothing : x)
    count(t) = (k = tryparse(Int, strip(shown_text(t))); k === nothing || k < 1 ? nothing : k)
    values() = (nanotube_indices(shown_text(indices)), num(bond), count(periods), num(vacuum),
                tube_elements(shown_text(elements)))
    update(_...) = (preview[] = nanotube_preview(values()...))
    foreach(t -> Makie.on(update, t.displayed_string), (indices, periods, bond, vacuum, elements))
    update()
    Makie.on(create.clicks) do _
        nm, b, k, vac, els = values()
        if any(isnothing, (nm, b, k, vac, els))
            v.status[] = nanotube_preview(nm, b, k, vac, els)
            return
        end
        c = nanotube(nm...; bond=b, periods=k, vacuum=vac, species=els)
        open_in_tab!(v, c)
        v.status[] = "Created the $(c.title): $(natoms(c)) atoms (not saved yet)"
    end
    finish_dialog!(v, :nanotube, "Nanotube", f)
    return (; indices, periods, bond, vacuum, elements, preview, create)
end

"""One element (both sublattices) or two, as symbols, or `nothing`."""
function tube_elements(s::AbstractString)
    t = split(replace(s, ',' => ' '))
    1 <= length(t) <= 2 || return nothing
    els = try
        [known_element(x) for x in t]
    catch
        return nothing
    end
    return length(els) == 1 ? (els[1], els[1]) : (els[1], els[2])
end

"""
    supercell_dialog(viewer) -> NamedTuple of the tool's controls

Supercell tool: repeat the cell n1 × n2 × n3, or give the new cell vectors
a′ b′ c′ as integer combinations of a, b, c (e.g. a′ = 1 1 0, b′ = -1 1 0 for
√2×√2 R45°). Presets fill in common rotated cells. A preview shows the atom
count and the new cell; Create replaces the structure (undoable).
"""
function supercell_dialog(v::Viewer)
    f = new_dialog(v, :supercell, "Supercell")
    f === nothing && return nothing
    heading!(f, "Supercell of the current cell")
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], "Type"; halign=:left, fontsize=PANEL_FONTSIZE)
    preset = Makie.Menu(f.g[r, 2]; options=first.(SUPERCELL_PRESETS), i_selected=1, fontsize=PANEL_FONTSIZE,
                        width=Makie.Relative(1))
    repeats = field!(f, "Repeats n1 n2 n3", "2 2 2", is_count3)
    heading!(f, "New cell vectors (in units of a, b, c)")
    vecs = [field!(f, label, "", is_int3) for label in ("a′ =", "b′ =", "c′ =")]
    preview = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], preview; halign=:left, fontsize=13, color=:gray25, justification=:left)
    note!(f, "Repeats fill in the vectors; edit them for any other cell.\n" *
             "Example: a′ = 1 1 0, b′ = -1 1 0, c′ = 0 0 1 is √2×√2 R45°.")
    create = button!(f, "Create supercell")

    # Repeats (and the preset) → the three vectors.
    function fill_vectors(_...)
        n = tryparse3(Int, shown_text(repeats))
        (n === nothing || any(<(1), n)) && return
        p = SUPERCELL_PRESETS[preset.i_selected[]][2]
        cols = p === nothing ? [(n[1], 0, 0), (0, n[2], 0), (0, 0, n[3])] :
               [n[1] .* p[1], n[2] .* p[2], (0, 0, n[3])]
        foreach((t, x) -> settext!(t, join(x, " ")), vecs, cols)
    end
    matrix() = supercell_matrix(shown_text.(vecs))
    update_preview(_...) = (preview[] = supercell_preview(v.crystal, matrix()))
    Makie.on(fill_vectors, repeats.displayed_string)
    Makie.on(fill_vectors, preset.i_selected)
    foreach(t -> Makie.on(update_preview, t.displayed_string), vecs)
    Makie.on(update_preview, v.header)                  # the structure changed
    fill_vectors(); update_preview()

    Makie.on(create.clicks) do _
        P = matrix()
        P === nothing && (v.status[] = "Enter a′ b′ c′ as three integers each"; return)
        what = isdiag(P) ? "make $(join(diag(P), "×")) supercell" :
               "make supercell a′ b′ c′ = " * join(("($(join(P[:, k], " ")))" for k in 1:3), " ")
        build!(c -> supercell(c, P), v, what)
    end
    finish_dialog!(v, :supercell, "Supercell", f)
    return (; preset, repeats, vecs, preview, create)
end

# ---------------------------------------------------------------------------

const SHIFT_COORDINATES = ["Fractional (of a, b, c)", "Cartesian (Å)"]
const SHIFT_ATOMS = ["All atoms", "Selected atoms"]

"""
Preview of a shift: how far, and how many atoms end up outside the cell
(or what is wrong).
"""
function shift_preview(c::Crystal, x, fractional::Bool, atoms, wrap_after::Bool)
    x === nothing && return "Enter the shift as three numbers"
    atoms !== nothing && isempty(atoms) && return "Select atoms first, or shift all atoms"
    n = atoms === nothing ? natoms(c) : length(atoms)
    f = fractional ? Vec3(x) : cart_to_frac(c.lattice, Vec3(x))
    r = frac_to_cart(c.lattice, f)
    head = @sprintf("Moves %s by %.3f Å\nΔx Δy Δz = %s Å    Δa Δb Δc = %s",
                    plural(n, "atom"), norm(r), vector_text(round.(r; digits=4)), vector_text(round.(f; digits=4)))
    wrap_after && return head * "   (then wrapped into the cell)"
    out = length(atoms_outside_cell(translate_atoms(c, f; fractional=true, atoms)))
    return head * "\n" * (out == 0 ? "All atoms stay inside the cell" :
                          "$(plural(out, "atom")) will be outside the cell")
end

"""
    shift_dialog(viewer) -> NamedTuple of the tool's controls

Shift atoms tool: translate all atoms (or the selected ones) by a vector in
fractional or Cartesian coordinates, optionally wrapping them into the cell
afterwards; "Selected atom → origin" fills in the shift that puts the selected
atom at (0, 0, 0). A second section wraps atoms that are outside the cell back
in by whole lattice vectors. Both are undoable.
"""
function shift_dialog(v::Viewer)
    f = new_dialog(v, :shift, "Shift atoms")
    f === nothing && return nothing
    heading!(f, "Shift atoms by a vector")
    menu_row(label, options, i) = (r = nextrow!(f);
        Makie.Label(f.g[r, 1], label; halign=:left, fontsize=PANEL_FONTSIZE);
        Makie.Menu(f.g[r, 2]; options, i_selected=i, fontsize=PANEL_FONTSIZE, width=Makie.Relative(1)))
    coords = menu_row("Coordinates", SHIFT_COORDINATES, 1)
    vector = field!(f, "Shift  x y z", "0 0 0", is_real3)
    atoms_menu = menu_row("Atoms", SHIFT_ATOMS, isempty(v.selection) ? 1 : 2)
    r = nextrow!(f)
    wrap_toggle = Makie.Toggle(f.g[r, 1]; active=false, halign=:left)
    Makie.Label(f.g[r, 2], "Wrap into the cell after shifting"; halign=:left, fontsize=PANEL_FONTSIZE)
    origin = button!(f, "Selected atom → origin")
    preview = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], preview; halign=:left, fontsize=13, color=:gray25, justification=:left)
    apply = button!(f, "Shift")

    heading!(f, "Wrap atoms into the cell")
    note!(f, "Moves atoms outside the cell back in by whole lattice vectors\n" *
             "(the structure and its bonds stay the same).")
    outside = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], outside; halign=:left, fontsize=13, color=:gray25)
    wrap_button = button!(f, "Wrap atoms into cell")

    fractional() = coords.i_selected[] == 1
    selected() = atoms_menu.i_selected[] == 2 ? selected_indices(v) : nothing
    function update(_...)
        c = v.crystal
        preview[] = shift_preview(c, tryparse3(Float64, shown_text(vector)), fractional(), selected(),
                                  wrap_toggle.active[])
        n = length(atoms_outside_cell(c; atoms=selected()))
        outside[] = (selected() === nothing ? "" : "Selected atoms: ") *
                    (n == 0 ? "none outside the cell" : "$(plural(n, "atom")) outside the cell")
    end
    # Switching coordinates converts the typed shift, so it stays the same shift.
    Makie.on(coords.i_selected) do i
        x = tryparse3(Float64, shown_text(vector))
        L = v.crystal.lattice
        x === nothing || settext!(vector, vector_text(i == 1 ? cart_to_frac(L, x) : frac_to_cart(L, x)))
        update()
    end
    Makie.on(update, vector.displayed_string)
    Makie.on(update, atoms_menu.i_selected)
    Makie.on(update, wrap_toggle.active)
    Makie.on(update, v.header)                  # the structure changed
    Makie.on(update, v.info)                    # the selection changed
    update()

    Makie.on(origin.clicks) do _
        idx = selected_indices(v)
        isempty(idx) && (v.status[] = "Select the atom to put at the origin"; return)
        a = first(selected_atoms(v))            # the copy that was clicked
        f0 = v.display.crystal.frac[a.index] + a.image
        settext!(vector, vector_text(fractional() ? -f0 : -frac_to_cart(v.crystal.lattice, f0)))
        atoms_menu.i_selected[] = 1
    end
    Makie.on(apply.clicks) do _
        x = tryparse3(Float64, shown_text(vector))
        x === nothing && (v.status[] = "Enter the shift as three numbers"; return)
        iszero(x) && (v.status[] = "Enter a non-zero shift"; return)
        atoms = selected()
        atoms !== nothing && isempty(atoms) && (v.status[] = "Select atoms first, or choose All atoms"; return)
        shift_atoms!(v, x; fractional=fractional(), atoms, wrap_after=wrap_toggle.active[])
    end
    Makie.on(_ -> wrap_atoms!(v; atoms=selected()), wrap_button.clicks)
    finish_dialog!(v, :shift, "Shift atoms", f)
    return (; coords, vector, atoms=atoms_menu, wrap=wrap_toggle, origin, preview, apply,
            outside, wrap_button)
end

# ---------------------------------------------------------------------------

const ROTATE_AXES = ["Cartesian  x y z", "Lattice  u v w", "Through 2 selected atoms", "Viewing direction"]
const ROTATE_CENTRES = ["Centre of the rotated atoms", "First selected atom", "Point (fractional)", "Point (Cartesian Å)"]

"""
The rotation chosen in the Rotate atoms tool, as a NamedTuple (angle,
Cartesian axis and centre, atoms and images to rotate, `about` text), or a
string saying what is missing. Selected atoms are taken as the copies that
were clicked, so a group across a cell face turns as one piece.
"""
function rotation_setup(v::Viewer; angle_text, axis_choice, axis_text, centre_choice, centre_text, selected_only)
    c, L = v.crystal, v.crystal.lattice
    angle = tryparse(Float64, strip(angle_text))
    angle === nothing && return "Enter the angle in degrees"
    sel = selected_atoms(v)                                  # drawn copies, in click order
    atoms = images = nothing
    if selected_only
        isempty(sel) && return "Select atoms first, or rotate all atoms"
        firsts = [a for (k, a) in enumerate(sel) if all(b -> b.index != a.index, sel[1:k-1])]
        atoms = [a.index for a in firsts]
        images = [Shift(round.(Int, cart_to_frac(L, a.position) - c.frac[a.index])) for a in firsts]
    end
    centre = nothing
    if axis_choice in (1, 2)
        x = tryparse3(Float64, axis_text)
        (x === nothing || iszero(x)) && return "Enter the axis as three numbers (not 0 0 0)"
        axis = axis_choice == 1 ? Vec3(x) : L.matrix * Vec3(x)
        about = axis_choice == 1 ? "[$(vector_text(x))]" : "[$(vector_text(x))] of the cell"
    elseif axis_choice == 3
        length(sel) >= 2 || return "Select two atoms for the axis (the first two you clicked)"
        axis, centre = sel[2].position - sel[1].position, sel[1].position
        about = "the $(display_label(c, sel[1]))–$(display_label(c, sel[2])) axis"
    else
        cam = camera(v.lscene)
        axis = Vec3(cam.eyeposition[] - cam.lookat[])
        about = "the viewing direction"
    end
    if centre === nothing
        if centre_choice == 1
            r = atoms === nothing ? cartesian_positions(c) :
                [frac_to_cart(L, c.frac[i] + s) for (i, s) in zip(atoms, images)]
            isempty(r) && return "Nothing to rotate: there are no atoms"
            centre = sum(r) / length(r)
        elseif centre_choice == 2
            isempty(sel) && return "Select the atom to rotate about"
            centre = sel[1].position
        else
            x = tryparse3(Float64, centre_text)
            x === nothing && return "Enter the centre as three numbers"
            centre = centre_choice == 3 ? frac_to_cart(L, Vec3(x)) : Vec3(x)
        end
    end
    return (; angle, axis, center=Vec3(centre), atoms, images, about)
end

"""Preview of a rotation: what turns about what, the closest contact afterwards, atoms outside the cell."""
function rotate_preview(c::Crystal, r, wrap_after::Bool)
    r isa AbstractString && return r
    n = r.atoms === nothing ? natoms(c) : length(r.atoms)
    lines = [@sprintf("Turns %s by %s° about %s\nthrough (%s) Å", plural(n, "atom"), vector_text([r.angle]),
                      r.about, vector_text(round.(r.center; digits=3)))]
    trial = try
        rotate_atoms(c, r.angle, r.axis; r.center, r.atoms, r.images)
    catch err
        return "Cannot rotate: " * sprint(showerror, err)
    end
    contact = closest_contact(trial, something(r.atoms, 1:natoms(c)))
    if contact !== nothing
        i, j, d = contact
        close = d < 0.6 * (element(c.species[i]).covalent_radius + element(c.species[j]).covalent_radius)
        push!(lines, @sprintf("Closest contact afterwards: %s–%s %.3f Å", atom_label(c, i), atom_label(c, j), d) *
                     (close ? "  (too close!)" : ""))
    end
    if wrap_after
        push!(lines, "Then wrapped into the cell")
    else
        out = length(atoms_outside_cell(trial))
        out > 0 && push!(lines, "$(plural(out, "atom")) will be outside the cell")
    end
    r.atoms === nothing && any(c.pbc) &&
        push!(lines, "The cell does not turn: unless this is a symmetry\noperation, the crystal changes.")
    return join(lines, "\n")
end

"""
    rotate_dialog(viewer) -> NamedTuple of the tool's controls

Rotate atoms tool: turn all atoms (or the selected ones) by an angle about
an axis (Cartesian, a lattice direction, the line through two selected atoms,
or the viewing direction) through a centre (the rotated atoms' centroid, the
first selected atom, or a point). The cell is not rotated. A preview names
the closest contact afterwards to catch overlaps. Undoable.
"""
function rotate_dialog(v::Viewer)
    f = new_dialog(v, :rotate, "Rotate atoms")
    f === nothing && return nothing
    heading!(f, "Rotate atoms about an axis")
    menu_row(label, options, i) = (r = nextrow!(f);
        Makie.Label(f.g[r, 1], label; halign=:left, fontsize=PANEL_FONTSIZE);
        Makie.Menu(f.g[r, 2]; options, i_selected=i, fontsize=PANEL_FONTSIZE, width=Makie.Relative(1)))
    atoms_menu = menu_row("Atoms", SHIFT_ATOMS, isempty(v.selection) ? 1 : 2)
    axis_menu = menu_row("Axis", ROTATE_AXES, 1)
    axis = field!(f, "Axis direction", "0 0 1", is_real3)
    centre_menu = menu_row("Centre", ROTATE_CENTRES, 1)
    centre = field!(f, "Centre point", "0.5 0.5 0.5", is_real3)
    angle = field!(f, "Angle (°)", "90", isnum)
    r = nextrow!(f)
    wrap_toggle = Makie.Toggle(f.g[r, 1]; active=false, halign=:left)
    Makie.Label(f.g[r, 2], "Wrap into the cell after rotating"; halign=:left, fontsize=PANEL_FONTSIZE)
    preview = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], preview; halign=:left, fontsize=13, color=:gray25, justification=:left)
    apply = button!(f, "Rotate")
    note!(f, "Positive angles turn counterclockwise when the axis points\n" *
             "toward you. Through 2 atoms: the line from the first atom you\n" *
             "clicked to the second; it also fixes the centre.")

    setup() = rotation_setup(v; angle_text=shown_text(angle), axis_choice=axis_menu.i_selected[],
                             axis_text=shown_text(axis), centre_choice=centre_menu.i_selected[],
                             centre_text=shown_text(centre), selected_only=atoms_menu.i_selected[] == 2)
    update(_...) = (preview[] = rotate_preview(v.crystal, setup(), wrap_toggle.active[]))
    foreach(t -> Makie.on(update, t.displayed_string), (axis, centre, angle))
    foreach(m -> Makie.on(update, m.i_selected), (atoms_menu, axis_menu, centre_menu))
    Makie.on(update, wrap_toggle.active)
    Makie.on(update, v.header)                  # the structure changed
    Makie.on(update, v.info)                    # the selection changed
    update()

    Makie.on(apply.clicks) do _
        r = setup()
        r isa AbstractString && (v.status[] = r; return)
        iszero(mod(r.angle, 360)) && (v.status[] = "Enter an angle that is not a multiple of 360°"; return)
        rotate_atoms!(v, r.angle, r.axis; r.center, r.atoms, r.images, r.about, wrap_after=wrap_toggle.active[])
    end
    finish_dialog!(v, :rotate, "Rotate atoms", f)
    return (; atoms=atoms_menu, axis_menu, axis, centre_menu, centre, angle, wrap=wrap_toggle,
            preview, apply)
end

# ---------------------------------------------------------------------------

const LIBRARY_KINDS = [("All structures", nothing), ("Crystals", :crystal), ("Molecules", :molecule)]

"""Group menu of the library tool: "All groups", then the groups of `entries` in library order."""
library_group_options(entries) = ["All groups"; unique(e.category for e in entries if !isempty(e.category))]

"""Whether a library entry matches the search text (name, formula, group or use; any case; all words)."""
function library_matches(e::LibraryEntry, query::AbstractString)
    text = lowercase(join((e.name, e.formula, e.category, e.use), " "))
    return all(w -> occursin(w, text), split(lowercase(query)))
end

"""Text about a library entry: what it is, its size or symmetry, and where it comes from."""
function library_info(e::LibraryEntry)
    c = try
        load_library(e)
    catch err
        return "Cannot read $(basename(e.file)): " * sprint(showerror, err)
    end
    shape = if e.kind === :crystal
        sg = try
            spacegroup(c)
        catch
            nothing
        end
        a, b, cc, α, β, γ = lattice_parameters(c.lattice)
        @sprintf("%s, %d atoms in the cell%s\na b c = %.4f %.4f %.4f Å   α β γ = %.1f %.1f %.1f°",
                 e.formula, natoms(c), sg === nothing ? "" : ", $(sg[2]) ($(sg[1]))", a, b, cc, α, β, γ)
    else
        "$(e.formula), $(plural(natoms(c), "atom")) (shown in a box with 5 Å around it)"
    end
    lines = [shape, isempty(e.category) ? "" : "Group: " * e.category, isempty(e.use) ? "" : "Used as: " * e.use,
             e.description, "Source: " * e.source, e.reference]
    isempty(e.doi) || push!(lines, "DOI: " * e.doi)
    isempty(e.note) || push!(lines, "Note: " * e.note)
    return join(filter(!isempty, lines), "\n")
end

"""
    open_library!(viewer, entry or name) -> Bool

Replace the structure with one from the library (undoable). It is Untitled
until saved, so Save never writes into the library. `box` (vacuum in Å) puts
a molecule into a periodic box for plane-wave codes (`cubic` box by default).
"""
open_library!(v::Viewer, name::AbstractString; kw...) = open_library!(v, library_entry(name); kw...)

function open_library!(v::Viewer, e::LibraryEntry; box::Union{Nothing,Real}=nothing, cubic::Bool=true)
    load(_) = (c = load_library(e); box === nothing || e.kind !== :molecule ? c : periodic_box(c; vacuum=box, cubic))
    what = "$(e.name) from the library" *
           (box !== nothing && e.kind === :molecule ? " in a $(vector_text([box])) Å vacuum box" : "")
    c = try
        load(nothing)
    catch err
        err isa InterruptException && rethrow()
        v.status[] = "Could not open $what: " * sprint(showerror, err)
        return false
    end
    open_in_tab!(v, c)
    v.status[] = "Opened $what: $(natoms(c)) atoms"
    return true
end

"""
    library_dialog(viewer) -> NamedTuple of the tool's controls

Structure library tool: choose crystals, molecules or all, narrow the list
by group (e.g. "Organic: hydrocarbons", "Inorganic: hydrides") or by a
search (name, formula or use such as "etch"), pick a structure and see
its formula, symmetry, cell and source (paper, database entry), then Open
it in the main window.
"""
function library_dialog(v::Viewer)
    width = 560
    f = new_dialog(v, :library, "Structure library"; width)
    f === nothing && return nothing
    heading!(f, "Open a known structure")
    note!(f, "Crystals: Crystallography Open Database (published structures).\n" *
             "Molecules: experimental gas-phase geometries from the literature.")
    menu_row(label, options) = (r = nextrow!(f);
        Makie.Label(f.g[r, 1], label; halign=:left, fontsize=PANEL_FONTSIZE);
        Makie.Menu(f.g[r, 2]; options, i_selected=1, fontsize=PANEL_FONTSIZE, width=Makie.Relative(1)))
    kind = menu_row("Show", first.(LIBRARY_KINDS))
    group = menu_row("Group", library_group_options(structure_library()))
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], "Search"; halign=:left, fontsize=PANEL_FONTSIZE)
    search = sync_display!(Makie.Textbox(f.g[r, 2]; placeholder="name, formula or use (e.g. etch)",
                                         width=Makie.Relative(1), fontsize=PANEL_FONTSIZE, textpadding=(8, 8, 5, 5)))
    push!(f.boxes, search)
    entries = Makie.Observable(structure_library())
    label(e) = "$(e.name)  ($(e.formula))"
    structure = menu_row("Structure", label.(entries[]))
    info = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], info; halign=:left, fontsize=13, color=:gray25, justification=:left)
    r = nextrow!(f)
    boxed = Makie.Toggle(f.g[r, 1]; active=false, halign=:left)
    Makie.Label(f.g[r, 2], "Molecules: in a periodic cubic box (for plane-wave codes)"; halign=:left,
                fontsize=PANEL_FONTSIZE)
    vacuum = field!(f, "Vacuum (Å)", "10", s -> isnum(s) && parse(Float64, s) > 0)
    open = button!(f, "Open")

    current() = (i = structure.i_selected[]; i === nothing || !(1 <= i <= length(entries[])) ? nothing : entries[][i])
    update(_...) = (e = current(); info[] = e === nothing ? "No structure matches" : library_info(e))
    query() = search.focused[] ? search.displayed_string[] : something(search.stored_string[], "")
    of_kind() = structure_library(; kind=LIBRARY_KINDS[something(kind.i_selected[], 1)][2])
    function refilter(_...)
        groups = group.options[]
        g = something(group.i_selected[], 1)
        chosen = 1 < g <= length(groups) ? groups[g] : nothing
        entries[] = filter(e -> (chosen === nothing || e.category == chosen) && library_matches(e, query()), of_kind())
        structure.options[] = isempty(entries[]) ? ["(none)"] : label.(entries[])
        structure.i_selected[] = 1
        update()
    end
    Makie.on(kind.i_selected) do _
        group.options[] = library_group_options(of_kind())
        group.i_selected[] = 1                                    # also refilters
    end
    Makie.on(refilter, group.i_selected)
    Makie.on(refilter, search.stored_string)
    Makie.on(s -> search.focused[] && refilter(), search.displayed_string)    # filter while typing
    Makie.on(update, structure.i_selected)
    update()
    Makie.on(open.clicks) do _
        e = current()
        e === nothing && (v.status[] = "Choose a structure first"; return)
        box = boxed.active[] ? tryparse(Float64, shown_text(vacuum)) : nothing
        boxed.active[] && (box === nothing || box <= 0) && (v.status[] = "Enter the vacuum in Å (a positive number)"; return)
        open_library!(v, e; box)
    end
    finish_dialog!(v, :library, "Structure library", f; labelwidth=0.22)
    return (; kind, group, search, structure, entries, info, boxed, vacuum, open)
end

# ---------------------------------------------------------------------------

"""Cell summary: periodicity, a b c α β γ, volume, density, vectors."""
function cell_info(c::Crystal)
    a, b, cc, α, β, γ = lattice_parameters(c.lattice)
    V = cell_volume(c.lattice)
    mass = sum((element(s).mass for s in c.species if element(s).Z > 0); init=0.0)
    kind = all(c.pbc) ? "Periodic crystal" : !any(c.pbc) ? "Molecule (not periodic): the box is only a container" :
           "Periodic along " * join(("a", "b", "c")[collect(c.pbc)], ", ")
    vecs = join((@sprintf("%s = (%.4f, %.4f, %.4f)", n, v...) for (n, v) in zip(("a", "b", "c"), cell_vectors(c.lattice))), "   ")
    return @sprintf("%s, %s\na b c = %.4f %.4f %.4f Å   α β γ = %.2f %.2f %.2f°\nV = %.3f Å³   density %.4f g/cm³\n%s",
                    kind, plural(natoms(c), "atom"), a, b, cc, α, β, γ, V, mass / V * 1.66053907, vecs)
end

"""Which exporters can express the cell constraints `cc` exactly."""
function cell_constraints_preview(c::Crystal, cc)
    cc = try
        check_cell_constraints(cc)
    catch err
        return sprint(showerror, err)
    end
    t = with_constraints(c; cell=cc)
    _, _, qe = cell_constraints_qe(t)
    _, vasp = cell_constraints_vasp(t)
    _, abinit = abinit_optcell(t)
    siesta = cc == (0, 0, 0, 0, 0, 0) ? "cell fixed" : cc == FREE_CELL ? "MD.VariableCell T" :
             "MD.VariableCell T relaxes the whole cell (constraints not applied)"
    dftk = cc == (0, 0, 0, 0, 0, 0) ? "variablecell=false" : cc == FREE_CELL ? "variablecell=true" :
           "variablecell=true relaxes the whole cell (constraints not applied)"
    return describe_cell_constraints(cc) * "\nQE: " * qe * "\nABINIT: " * abinit * "\nSIESTA: " * siesta *
           "\nDFTK: " * dftk *
           "\nVASP: " * vasp * "\nCASTEP: CELL_CONSTRAINTS " * join(cc[1:3], " ") * " / " * join(cc[4:6], " ")
end

"""Atom constraint summary: how many atoms are fixed, and the selection's state."""
function atom_constraints_info(v::Viewer)
    c = v.crystal
    full, part = count(all, c.fixed), count(f -> any(f) && !all(f), c.fixed)
    s = full + part == 0 ? "No atom is fixed" :
        "Fixed: " * plural(full, "atom") * (part > 0 ? ", partly fixed: " * plural(part, "atom") : "")
    idx = selected_indices(v)
    isempty(idx) && return s * "\nSelect atoms in the main window to fix them"
    nf = count(i -> any(c.fixed[i]), idx)
    return s * "\nSelected: " * plural(length(idx), "atom") * (nf > 0 ? ", $nf of them fixed" : ", none fixed")
end

"""
    cell_dialog(viewer) -> NamedTuple of the tool's controls

Cell and constraints tool:
- Cell: periodicity, lattice parameters, volume, density and vectors.
- Periodic box: put a molecule in a box with a given vacuum (for plane-wave
  codes), or treat a structure as a molecule.
- Lattice constraints for geometry optimization: presets or numbers for
  a b c α β γ (0 fixed, equal numbers tied) with how each code gets them
  (Quantum ESPRESSO, ABINIT, SIESTA, DFTK, VASP, CASTEP).
- Atom constraints: fix the fractional x, y, z of the selected atoms.
All changes are undoable and are written by the exporters.
"""
function cell_dialog(v::Viewer)
    width = 560
    f = new_dialog(v, :cell, "Cell and constraints"; width)
    f === nothing && return nothing
    small_label(obs) = Makie.Label(f.g[nextrow!(f), 1:2], obs; halign=:left, fontsize=13, color=:gray25,
                                   justification=:left)
    heading!(f, "Cell")
    info = Makie.Observable("")
    small_label(info)

    heading!(f, "Periodic box (molecules)")
    vacuum = field!(f, "Vacuum (Å)", "10", s -> isnum(s) && parse(Float64, s) > 0)
    r = nextrow!(f)
    cubic = Makie.Toggle(f.g[r, 1]; active=true, halign=:left)
    Makie.Label(f.g[r, 2], "Cubic box"; halign=:left, fontsize=PANEL_FONTSIZE)
    box_label = Makie.Observable("Put the molecule in a periodic box")
    box = button!(f, box_label)

    heading!(f, "Lattice constraints (geometry optimization)")
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], "Preset"; halign=:left, fontsize=PANEL_FONTSIZE)
    preset = Makie.Menu(f.g[r, 2]; options=[first.(CELL_CONSTRAINT_PRESETS); "Custom"], i_selected=1,
                        fontsize=PANEL_FONTSIZE, width=Makie.Relative(1))
    lengths = field!(f, "a b c", "1 2 3", is_count3)
    angles = field!(f, "α β γ", "4 5 6", is_count3)
    note!(f, "0 = fixed, equal numbers = vary together, different = independent\n(CASTEP CELL_CONSTRAINTS).")
    cell_preview = Makie.Observable("")
    small_label(cell_preview)
    apply_cell = button!(f, "Apply lattice constraints")

    heading!(f, "Atom constraints (fixed fractional coordinates)")
    atoms_info = Makie.Observable("")
    small_label(atoms_info)
    tg = Makie.GridLayout(f.g[nextrow!(f), 1:2]; halign=:left)
    axis_toggles = map(1:3) do k
        t = Makie.Toggle(tg[1, 2k-1]; active=true)
        Makie.Label(tg[1, 2k], "Fix $(("x", "y", "z")[k])"; fontsize=PANEL_FONTSIZE, padding=(4, 20, 0, 0))
        t
    end
    bg = Makie.GridLayout(f.g[nextrow!(f), 1:2])
    mk(i, label) = with_click_feedback!(Makie.Button(bg[1, i]; label, width=Makie.Relative(1), fontsize=PANEL_FONTSIZE,
                                                     padding=(6, 6, 5, 5)))
    fix, free, free_all = mk(1, "Fix selected"), mk(2, "Free selected"), mk(3, "Free all")
    note!(f, "x y z are fractional coordinates along a, b, c. Exported as VASP selective\n" *
             "dynamics, CASTEP IONIC_CONSTRAINTS and QE if_pos. Fixed atoms show × (+ partly).")

    typed_cc() = (x = tryparse3(Int, shown_text(lengths)); y = tryparse3(Int, shown_text(angles));
                  x === nothing || y === nothing ? nothing : (x..., y...))
    function show_cc(cc)
        settext!(lengths, join(cc[1:3], " ")); settext!(angles, join(cc[4:6], " "))
        k = findfirst(p -> p[2] == cc, CELL_CONSTRAINT_PRESETS)
        preset.i_selected[] == something(k, length(CELL_CONSTRAINT_PRESETS) + 1) ||
            (preset.i_selected[] = something(k, length(CELL_CONSTRAINT_PRESETS) + 1))
    end
    function update(_...)
        c = v.crystal
        info[] = cell_info(c)
        box_label[] = any(c.pbc) ? "Make it a molecule (not periodic)" : "Put the molecule in a periodic box"
        cc = typed_cc()
        cell_preview[] = cc === nothing ? "Enter three numbers ≥ 0 for a b c and for α β γ" :
                         cell_constraints_preview(c, cc)
        # The preset menu follows typed numbers (Custom when none matches).
        if cc !== nothing
            k = something(findfirst(p -> p[2] == Tuple(cc), CELL_CONSTRAINT_PRESETS), length(CELL_CONSTRAINT_PRESETS) + 1)
            preset.i_selected[] == k || (preset.i_selected[] = k)
        end
        atoms_info[] = atom_constraints_info(v)
    end
    Makie.on(preset.i_selected) do i
        i <= length(CELL_CONSTRAINT_PRESETS) && show_cc(CELL_CONSTRAINT_PRESETS[i][2])
    end
    foreach(t -> Makie.on(update, t.displayed_string), (lengths, angles))
    Makie.on(_ -> (show_cc(v.crystal.cell_constraints); update()), v.header)     # structure changed (also undo)
    Makie.on(update, v.info)                                                         # selection changed
    show_cc(v.crystal.cell_constraints)
    update()

    Makie.on(box.clicks) do _
        if any(v.crystal.pbc)
            make_molecule!(v)
        else
            x = tryparse(Float64, shown_text(vacuum))
            (x === nothing || x <= 0) && (v.status[] = "Enter the vacuum in Å (a positive number)"; return)
            periodic_box!(v; vacuum=x, cubic=cubic.active[])
        end
    end
    Makie.on(apply_cell.clicks) do _
        cc = typed_cc()
        cc === nothing ? (v.status[] = "Enter three numbers ≥ 0 for a b c and for α β γ") : set_cell_constraints!(v, cc)
    end
    Makie.on(_ -> fix_selected!(v, Tuple(t.active[] for t in axis_toggles)), fix.clicks)
    Makie.on(_ -> free_selected!(v), free.clicks)
    Makie.on(_ -> free_all_atoms!(v), free_all.clicks)
    finish_dialog!(v, :cell, "Cell and constraints", f; labelwidth=0.25)
    return (; info, vacuum, cubic, box, preset, lengths, angles, cell_preview, apply_cell, atoms_info,
            axis_toggles, fix, free, free_all)
end

# ---------------------------------------------------------------------------

const BRIDGE_OFFSET = 0.5     # Å, default offset of a bridge from the bond midpoint

"""Element, offset and side for a bridge: from the Oxygen bridge tool when it is open, else defaults."""
function bridge_settings(v::Viewer)
    d = get(v.dialogs, :bridge_form, nothing)
    d === nothing && return (; element=:O, offset=BRIDGE_OFFSET, flip=false)
    el = shown_text(d.element)
    x = tryparse(Float64, shown_text(d.offset))
    return (; element=is_element(el) ? known_element(el) : :O, offset=something(x, BRIDGE_OFFSET), flip=d.flip.active[])
end

"""
Preview of a bridge for the current selection: the bond, where the atom
goes, its distances to the two bond atoms and the bridge angle, the closest
other atom, and whether it is inside the triangle.
"""
function bridge_preview(v::Viewer, element, offset, flip)
    sel = bridge_atoms(v)
    sel isa AbstractString && return sel
    offset === nothing && return "Enter the offset in Å (a number ≥ 0)"
    element === nothing && return "Enter the element of the bridging atom (e.g. O)"
    dc = v.display.crystal
    a, b, third = (x.position for x in sel)
    na, nb, nc = (display_label(dc, x) for x in sel)
    p = try
        bridge_position(a, b, third; offset, flip)
    catch err
        return "Cannot place a bridge: " * sprint(showerror, err)
    end
    lines = [@sprintf("Bond %s–%s (%.3f Å); plane and side from %s", na, nb, norm(b - a), nc),
             @sprintf("%s at (%.3f, %.3f, %.3f) Å: %s–%s %.3f Å, %s–%s %.3f Å, %s–%s–%s %.1f°",
                      element, p..., na, element, norm(p - a), nb, element, norm(p - b), na, element, nb,
                      bond_angle(a, p, b))]
    trial = try
        add_atom(v.crystal, element, p; cartesian=true, min_distance=0)
    catch err
        return join([lines; "Cannot place a bridge: " * sprint(showerror, err)], "\n")
    end
    n = natoms(trial)
    others = [j for j in 1:n-1 if j ∉ (sel[1].index, sel[2].index)]
    if !isempty(others)
        d, j = findmin(j -> min_image_distance(trial, n, j), others)
        push!(lines, d < 0.1 ? "$(atom_label(trial, others[j])) is already there (Flip last bridge moves it)" :
                     @sprintf("Closest other atom: %s at %.3f Å", atom_label(trial, others[j]), d) *
                     (d < 0.5 ? "  (too close: cannot add)" : ""))
    end
    limit = bridge_max_offset(a, b, third)
    push!(lines, flip ? "On the other side of the bond (−v), away from $nc" :
                 offset <= limit ? @sprintf("Inside the triangle (offsets up to %.2f Å stay inside)", limit) :
                 @sprintf("Outside the triangle: the offset is larger than %.2f Å", limit))
    return join(lines, "\n")
end

"""
    bridge_dialog(viewer) -> NamedTuple of the tool's controls

Oxygen bridge tool: select three atoms in the main window (the bond's two
atoms, then a third). The bridging atom (O by default) goes at the bond's
middle plus an offset v perpendicular to the bond, in the plane of the three
atoms, towards the third (inside the triangle), or with "Other side" at −v.
"Flip last bridge" moves the bridge just added to the other side. Undoable.
"""
function bridge_dialog(v::Viewer)
    width = 520
    f = new_dialog(v, :bridge_tool, "Oxygen bridge"; width)
    f === nothing && return nothing
    heading!(f, "Bridge a bond with oxygen")
    note!(f, "Select three atoms (Shift+click): the first two are the bond, the third\n" *
             "sets the plane and the side. The atom goes to the bond's middle, moved by\n" *
             "the offset towards the third atom (inside the triangle).")
    element = field!(f, "Element", "O", is_element)
    offset = field!(f, "Offset from the bond (Å)", string(BRIDGE_OFFSET), s -> isnum(s) && parse(Float64, s) >= 0)
    r = nextrow!(f)
    flip = Makie.Toggle(f.g[r, 1]; active=false, halign=:left)
    Makie.Label(f.g[r, 2], "Other side of the bond (−v)"; halign=:left, fontsize=PANEL_FONTSIZE)
    preview = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], preview; halign=:left, fontsize=13, color=:gray25, justification=:left)
    bg = Makie.GridLayout(f.g[nextrow!(f), 1:2])
    mk(i, label) = with_click_feedback!(Makie.Button(bg[1, i]; label, width=Makie.Relative(1), fontsize=PANEL_FONTSIZE,
                                                     padding=(6, 6, 5, 5)))
    add, flip_last = mk(1, "Add bridge"), mk(2, "Flip last bridge")

    form = (; element, offset, flip)
    v.dialogs[:bridge_form] = form
    function update(_...)
        s = bridge_settings(v)
        el, x = shown_text(element), tryparse(Float64, shown_text(offset))
        preview[] = bridge_preview(v, is_element(el) ? known_element(el) : nothing,
                                   x === nothing || x < 0 ? nothing : x, s.flip)
    end
    foreach(t -> Makie.on(update, t.displayed_string), (element, offset))
    Makie.on(update, flip.active)
    Makie.on(update, v.info)                    # the selection changed
    Makie.on(update, v.header)                  # the structure changed
    update()

    Makie.on(add.clicks) do _
        el, x = shown_text(element), tryparse(Float64, shown_text(offset))
        is_element(el) || (v.status[] = "Enter the element of the bridging atom (e.g. O)"; return)
        (x === nothing || x < 0) && (v.status[] = "Enter the offset in Å (a number ≥ 0)"; return)
        add_bridge!(v; element=known_element(el), offset=x, flip=flip.active[])
    end
    Makie.on(_ -> flip_bridge!(v), flip_last.clicks)
    finish_dialog!(v, :bridge_tool, "Oxygen bridge", f; labelwidth=0.42)
    return (; element, offset, flip, preview, add, flip_last)
end
