# Open and Save as tools: a small file browser drawn with Makie (WSL has no
# native file dialog). Folders and readable structure files are listed; click
# selects and previews, double-click (or Open / Enter) opens a folder or loads
# a file into the main window as an undoable edit. In the Save as window a
# click on a file copies its name into the File name box, and Save writes the
# structure (or an image of the view) in the chosen format.

const BROWSER_ROWS = 16
const BROWSER_WIDTH = 520
const BROWSER_NAME_CHARS = 52
const DOUBLE_CLICK = 0.5       # s

const FORMAT_LABELS = Dict(:cif => "CIF", :xyz => "XYZ", :poscar => "POSCAR", :xdatcar => "VASP XDATCAR",
                           :vasprun => "vasprun.xml", :castep_cell => "CASTEP cell", :castep_geom => "CASTEP geom",
                           :castep_out => "CASTEP output",
                           :qe_in => "QE input", :qe_out => "QE output", :aims => "FHI-aims",
                           :siesta => "SIESTA fdf", :abinit => "ABINIT", :xsf => "XSF", :cube => "cube",
                           :pdb => "PDB", :lammps => "LAMMPS data", :car => "Materials Studio car",
                           :xsd => "Materials Studio xsd", :res => "SHELX res")

const EXAMPLES_DIR = normpath(joinpath(@__DIR__, "..", "..", "examples", "structures"))

"""
    browser_entries(dir; all_files=false) -> Vector of (name, isdir)

Folders first, then files, each sorted case-insensitively. Hidden entries are
skipped; files are limited to formats [`load_structure`](@ref) reads unless
`all_files`.
"""
function browser_entries(dir::AbstractString; all_files::Bool=false)
    entries = Tuple{String,Bool}[]
    for name in readdir(dir)
        startswith(name, ".") && continue
        isd = try isdir(joinpath(dir, name)) catch; false end
        (isd || all_files || structure_format(name) !== nothing) && push!(entries, (name, isd))
    end
    return sort!(entries; by=e -> (!e[2], lowercase(e[1])))
end

shorten(s, n) = length(s) <= n ? s : first(s, n - 1) * "…"

mutable struct FileBrowser
    v::Viewer
    figure::Makie.Figure
    dir::String
    entries::Vector{Tuple{String,Bool}}
    offset::Int                                   # entries scrolled past
    selected::Int                                 # index into entries, 0 = none
    all_files::Bool
    names::Makie.Observable{Vector{String}}       # per visible row
    kinds::Makie.Observable{Vector{String}}
    fills::Makie.Observable{Vector{RGBAf}}
    list::Any                                     # Axis the rows are drawn in (row k = y in [k-1, k])
    folder::Any                                   # folder text box
    position::Makie.Observable{String}            # "1–16 of 40"
    preview::Makie.Observable{String}             # about the selected entry
    last_click::Tuple{Int,Float64}
    mode::Symbol                                  # :open or :save
    name::Any                                     # File name text box (save mode)
    format::Any                                   # format Menu (save mode)
    confirm::String                               # existing file the next Save will replace
end

function refresh!(b::FileBrowser)
    n = length(b.entries)
    names, kinds, fills = String[], String[], RGBAf[]
    for k in 1:BROWSER_ROWS
        i = b.offset + k
        if i <= n
            name, isd = b.entries[i]
            push!(names, shorten(isd ? name * "/" : name, BROWSER_NAME_CHARS))
            push!(kinds, isd ? "folder" : get(FORMAT_LABELS, structure_format(name), ""))
        else
            push!(names, ""); push!(kinds, "")
        end
        push!(fills, i == b.selected ? RGBAf(0.26, 0.52, 0.96, 0.28) :
                     isodd(k) || i > n ? RGBAf(0, 0, 0, 0) : RGBAf(0, 0, 0, 0.035))
    end
    b.names[], b.kinds[], b.fills[] = names, kinds, fills
    b.position[] = n == 0 ? (b.all_files ? "Empty folder" : "No folders or structure files here") :
                   "$(b.offset + 1)–$(min(n, b.offset + BROWSER_ROWS)) of $n"
    return b
end

function scroll!(b::FileBrowser, by::Integer)
    b.offset = clamp(b.offset + by, 0, max(0, length(b.entries) - BROWSER_ROWS))
    refresh!(b)
end

"""Show folder `dir`; returns false (with a message) if it can't be listed."""
function navigate!(b::FileBrowser, dir::AbstractString)
    dir = normpath(abspath(expanduser(dir)))
    length(dir) > 1 && (dir = rstrip(dir, '/'))     # "a/b/" and "a/b" are the same folder
    entries = try
        browser_entries(dir; b.all_files)
    catch err
        err isa InterruptException && rethrow()
        b.v.status[] = "Cannot open folder $dir: " * sprint(showerror, err)
        return false
    end
    b.dir, b.entries, b.offset, b.selected = dir, entries, 0, 0
    b.folder.displayed_string[] = dir
    b.folder.stored_string[] != dir && (b.folder.stored_string[] = dir)
    b.preview[] = ""
    refresh!(b)
    return true
end

"""Select entry `i` (scrolling it into view) and describe it."""
function select_entry!(b::FileBrowser, i::Integer)
    1 <= i <= length(b.entries) || return
    b.selected = i
    i <= b.offset && (b.offset = i - 1)
    i > b.offset + BROWSER_ROWS && (b.offset = i - BROWSER_ROWS)
    refresh!(b)
    name, isd = b.entries[i]
    b.mode === :save && !isd && settext!(b.name, name)
    b.preview[] = isd ? "Folder — double-click or Open to enter" : try
        c = load_structure(joinpath(b.dir, name))
        a, bb, cc, α, β, γ = lattice_parameters(c.lattice)
        @sprintf("%s, %d atoms   a b c = %.3f %.3f %.3f Å   α β γ = %.1f %.1f %.1f°",
                 formula(c), natoms(c), a, bb, cc, α, β, γ)
    catch err
        err isa InterruptException && rethrow()
        "Cannot read: " * shorten(sprint(showerror, err), 120)
    end
end

"""Enter the selected folder or load the selected file into the viewer."""
function activate!(b::FileBrowser, i::Integer=b.selected)
    b.mode === :save && i == 0 && return save_from_browser!(b)
    1 <= i <= length(b.entries) || (b.v.status[] = "Select a file or folder first"; return false)
    name, isd = b.entries[i]
    path = joinpath(b.dir, name)
    isd && return navigate!(b, path)
    b.mode === :save && return save_from_browser!(b)
    ok = open_structure!(b.v, path)
    ok && (b.v.folder = b.dir)
    return ok
end

"""
    open_structure!(viewer, path) -> Bool

Open a structure file in a new tab (or in place of an untouched empty
cell; a file already open in a tab is shown there).
"""
function open_structure!(v::Viewer, path::AbstractString)
    c = try
        load_structure(path)
    catch err
        err isa InterruptException && rethrow()
        v.status[] = "Could not open $(basename(path)): " * sprint(showerror, err)
        return false
    end
    open_in_tab!(v, c; path)
    v.status[] = "Opened $(basename(path)): $(natoms(c)) atoms"
    return true
end

"""Save as: write to the folder shown and the typed file name in the chosen format."""
function save_from_browser!(b::FileBrowser)
    v = b.v
    name = shown_text(b.name)
    isempty(name) && (v.status[] = "Type a file name first"; return false)
    occursin('/', name) && (v.status[] = "Type a file name without a folder (choose the folder above)"; return false)
    path = joinpath(b.dir, with_format_extension(name, SAVE_FORMATS[b.format.i_selected[]][2]))
    if isfile(path) && b.confirm != path
        b.confirm = path
        v.status[] = "Replace $(basename(path))? It already exists: click Save again to replace it"
        return false
    end
    b.confirm = ""
    ok = save_to!(v, path)
    if ok
        v.folder = b.dir
        navigate!(b, b.dir)                   # show the new file
    end
    return ok
end

"""Visible row (1 = top) under the mouse, or 0 (outside the list or past the last entry)."""
function row_at(b::FileBrowser)
    Makie.is_mouseinside(b.list.scene) || return 0
    y = Makie.mouseposition(b.list.scene)[2]
    k = clamp(floor(Int, y) + 1, 1, BROWSER_ROWS)
    return b.offset + k <= length(b.entries) ? k : 0
end

"""Pixel position (figure coordinates) of the middle of visible row `k`, e.g. for simulated clicks."""
function row_pixel(b::FileBrowser, k::Integer)
    r = b.list.scene.viewport[]
    return (r.origin[1] + r.widths[1] / 2, r.origin[2] + r.widths[2] * (1 - (k - 0.5) / BROWSER_ROWS))
end

function install_browser_events!(b::FileBrowser)
    ev = Makie.events(b.figure.scene)
    Makie.on(ev.mousebutton; priority=100) do e
        (e.button == Makie.Mouse.left && e.action == Makie.Mouse.press) || return Makie.Consume(false)
        k = row_at(b)
        k == 0 && return Makie.Consume(false)
        i = b.offset + k
        # Timed from the end of the previous click's handling, so a slow
        # preview (first load compiles the reader) doesn't eat the double-click.
        if b.last_click[1] == i && time() - b.last_click[2] < DOUBLE_CLICK
            b.last_click = (0, 0.0)
            activate!(b, i)
        else
            select_entry!(b, i)
            b.last_click = (i, time())
        end
        return Makie.Consume(true)
    end
    Makie.on(ev.scroll; priority=100) do (_, dy)
        Makie.is_mouseinside(b.list.scene) || return Makie.Consume(false)
        scroll!(b, dy > 0 ? -3 : 3)
        return Makie.Consume(true)
    end
    Makie.on(ev.keyboardbutton) do e
        e.action in (Makie.Keyboard.press, Makie.Keyboard.repeat) || return Makie.Consume(false)
        (b.folder.focused[] || (b.name !== nothing && b.name.focused[])) && return Makie.Consume(false)
        K = Makie.Keyboard
        n = length(b.entries)
        if e.key == K.down
            select_entry!(b, min(n, b.selected + 1))
        elseif e.key == K.up
            select_entry!(b, max(1, b.selected - 1))
        elseif e.key == K.page_down
            select_entry!(b, min(n, b.selected + BROWSER_ROWS))
        elseif e.key == K.page_up
            select_entry!(b, max(1, b.selected - BROWSER_ROWS))
        elseif e.key in (K.enter, K.kp_enter)
            activate!(b)
        elseif e.key == K.backspace
            navigate!(b, dirname(b.dir))
        else
            return Makie.Consume(false)
        end
        return Makie.Consume(true)
    end
end

"""
    open_dialog(viewer; show=true, dir=viewer.folder) -> FileBrowser

Open tool for loading structures into `viewer`. Lists folders and files
that [`load_structure`](@ref) reads (CIF, POSCAR/CONTCAR/.vasp, XYZ, CASTEP
.cell/.geom); "All files" lists everything. Keys: ↑/↓, Page Up/Down, Enter
opens, Backspace goes to the parent folder.
"""
open_dialog(v::Viewer; show::Bool=true, dir::AbstractString=v.folder) = file_dialog(v, :open; show, dir)

"""
    save_dialog(viewer; show=true, dir, format=nothing) -> FileBrowser

Save as tool: choose a folder, type a file name and pick a format (POSCAR,
CIF, XYZ, CASTEP .cell, Quantum ESPRESSO .in, or a PNG image of the view).
Replacing an existing file needs a second click on Save.
"""
function save_dialog(v::Viewer; show::Bool=true, dir::AbstractString=default_save_folder(v),
                     format::Union{Nothing,AbstractString}=nothing)
    b = file_dialog(v, :save; show, dir)
    b === nothing && return nothing
    ext = format !== nothing ? format : v.path === nothing ? ".vasp" : save_extension(v.path)
    k = something(findfirst(f -> f[2] == ext, SAVE_FORMATS), 1)
    b.format.i_selected[] = k
    stem = v.path === nothing ? (natoms(v.crystal) == 0 ? "structure" : formula(v.crystal)) :
           splitext(basename(v.path))[1]
    settext!(b.name, with_format_extension(stem, SAVE_FORMATS[k][2]))
    return b
end

function file_dialog(v::Viewer, mode::Symbol; show::Bool=true, dir::AbstractString=v.folder)
    saving = mode === :save
    title = saving ? "Save as" : "Open structure"
    d = new_dialog(v, mode, title, saving ? 760 : 690; show, width=BROWSER_WIDTH)
    d === nothing && return nothing
    fig, f = d

    heading!(f, title)
    folder = field!(f, "Folder", dir, s -> isdir(expanduser(strip(s))))
    r = nextrow!(f)
    nav = Makie.GridLayout(f.g[r, 1:2]; colgap=6)
    navbtn(k, label) = with_click_feedback!(Makie.Button(nav[1, k]; label, width=Makie.Relative(1),
                                                         fontsize=PANEL_FONTSIZE, padding=(6, 6, 4, 4)))
    up, home, here, examples = navbtn(1, "Up"), navbtn(2, "Home"), navbtn(3, "Working folder"), navbtn(4, "Examples")

    # The list: an Axis without decorations, one unit of y per row, top row first.
    list = Makie.Axis(f.g[nextrow!(f), 1:2]; height=22 * BROWSER_ROWS, yreversed=true,
                      limits=(0, 1, 0, BROWSER_ROWS), backgroundcolor=:transparent)
    Makie.hidedecorations!(list)
    Makie.hidespines!(list)
    foreach(i -> Makie.deregister_interaction!(list, i), (:rectanglezoom, :dragpan, :scrollzoom, :limitreset))
    names = Makie.Observable(fill("", BROWSER_ROWS))
    kinds = Makie.Observable(fill("", BROWSER_ROWS))
    fills = Makie.Observable(fill(RGBAf(0, 0, 0, 0), BROWSER_ROWS))
    Makie.poly!(list, [Makie.Rect2f(0, k - 1, 1, 1) for k in 1:BROWSER_ROWS]; color=fills, strokewidth=0)
    Makie.text!(list, [Makie.Point2f(0.015, k - 0.5) for k in 1:BROWSER_ROWS]; text=names,
                align=(:left, :center), fontsize=13)
    Makie.text!(list, [Makie.Point2f(0.985, k - 0.5) for k in 1:BROWSER_ROWS]; text=kinds,
                align=(:right, :center), fontsize=12, color=:gray45)

    r = nextrow!(f)
    position = Makie.Observable("")
    Makie.Label(f.g[r, 1], position; halign=:left, fontsize=12, color=:gray40)
    tg = Makie.GridLayout(f.g[r, 2]; colgap=6, halign=:right)
    all_toggle = Makie.Toggle(tg[1, 1]; active=false)
    Makie.Label(tg[1, 2], "All files"; fontsize=PANEL_FONTSIZE)
    preview = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], preview; halign=:left, fontsize=12, color=:gray25,
                word_wrap=true, width=BROWSER_WIDTH - 40)
    name = format = nothing
    if saving
        name = field!(f, "File name", "", s -> !isempty(strip(s)))
        r = nextrow!(f)
        Makie.Label(f.g[r, 1], "Format"; halign=:left, fontsize=PANEL_FONTSIZE)
        format = Makie.Menu(f.g[r, 2]; options=first.(SAVE_FORMATS), i_selected=1, fontsize=PANEL_FONTSIZE,
                            width=Makie.Relative(1))
    end
    go = button!(f, saving ? "Save" : "Open")

    b = FileBrowser(v, fig, "", Tuple{String,Bool}[], 0, 0, false, names, kinds, fills, list, folder,
                    position, preview, (0, 0.0), mode, name, format, "")
    if saving
        # The extension follows the format.
        Makie.on(format.i_selected) do k
            s = shown_text(name)
            isempty(s) || settext!(name, with_format_extension(s, SAVE_FORMATS[k][2]))
        end
        Makie.on(_ -> (b.confirm = ""), name.displayed_string)
        Makie.on(_ -> save_from_browser!(b), go.clicks)
    else
        Makie.on(_ -> activate!(b), go.clicks)
    end
    Makie.on(s -> s !== nothing && rstrip(normpath(abspath(expanduser(strip(s)))), '/') != b.dir &&
                  navigate!(b, strip(s)), folder.stored_string)
    Makie.on(_ -> navigate!(b, dirname(b.dir)), up.clicks)
    Makie.on(_ -> navigate!(b, homedir()), home.clicks)
    Makie.on(_ -> navigate!(b, pwd()), here.clicks)
    Makie.on(_ -> navigate!(b, EXAMPLES_DIR), examples.clicks)
    Makie.on(all_toggle.active) do x
        b.all_files = x
        navigate!(b, b.dir)
    end
    install_browser_events!(b)
    navigate!(b, isdir(dir) ? dir : pwd())

    finish_dialog!(v, mode, title, fig, f; show, labelwidth=saving ? 0.2 : 0.15)
    return b
end
