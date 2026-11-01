# Menu bar of the main window (File / Edit / Build / View), saving, recent
# files and the New structure tool.
#
# The current file is shown at the right of the menu bar, with "modified"
# after unsaved changes. Save writes back to that file when it is a POSCAR or
# XYZ file, or one AtomBond saved; a CIF or CASTEP file opened from disk is
# not overwritten without asking (Save opens Save as), since AtomBond writes
# only the cell and atoms and would drop the rest of such a file.

const SAVE_FORMATS = [("POSCAR (.vasp)", ".vasp"), ("CIF (.cif)", ".cif"), ("XYZ (.xyz)", ".xyz"),
                      ("CASTEP cell (.cell)", ".cell"), ("Quantum ESPRESSO pw.x (.in)", ".in"),
                      ("FHI-aims (geometry.in)", "geometry.in"), ("XSF, XCrySDen / VESTA (.xsf)", ".xsf"),
                      ("PDB (.pdb)", ".pdb"), ("LAMMPS data (.lmp)", ".lmp"), ("Materials Studio (.car)", ".car"),
                      ("SHELX / AIRSS (.res)", ".res"), ("PNG image of the view (.png)", ".png")]
const SAVE_EXTENSIONS = last.(SAVE_FORMATS)
const MAX_RECENT = 8
const MENUBAR_MENU_WIDTH = 290

# Files AtomBond wrote in this session: Save may overwrite them.
const WRITTEN = Set{String}()

"""Settings folder: \$ATOMBOND_HOME, or ~/.atombond."""
atombond_home() = get(ENV, "ATOMBOND_HOME", joinpath(homedir(), ".atombond"))
recent_list_file() = joinpath(atombond_home(), "recent.txt")

"""Recently opened or saved structure files that still exist, newest first."""
function recent_files()
    p = recent_list_file()
    isfile(p) || return String[]
    return filter(isfile, readlines(p))
end

function add_recent!(path::AbstractString)
    structure_format(path) === nothing && return   # only files that can be opened again
    path = abspath(path)
    list = [path; filter(!=(path), recent_files())]
    try
        mkpath(atombond_home())
        write(recent_list_file(), join(list[1:min(end, MAX_RECENT)], "\n") * "\n")
    catch err
        err isa InterruptException && rethrow()     # the list is a convenience only
    end
    return
end

"""File name for the menu bar: "GaAs.vasp — modified", "Untitled — not saved"."""
function update_document!(v::Viewer)
    v.document[] = v.path === nothing ? "Untitled — not saved" :
                   basename(v.path) * (get(v.saved, v.path, nothing) === v.crystal ? "" : " — modified")
    refresh_tab_labels!(v)
end

"""`name` with the extension of a save format (replacing a known one; POSCAR/CONTCAR names stay)."""
function with_format_extension(name::AbstractString, ext::AbstractString)
    ext == "geometry.in" && return "geometry.in"          # FHI-aims reads this name only
    stem, e = splitext(name)
    ext == ".vasp" && isempty(e) && occursin(r"^(POSCAR|CONTCAR)"i, name) && return String(name)
    lowercase(e) in SAVE_EXTENSIONS || lowercase(e) == ".poscar" || (stem = name)
    return stem * ext
end

# Format to save a file read in format `fmt` as (outputs as the matching input).
const SAVE_EXTENSION_OF = Dict(:cif => ".cif", :xyz => ".xyz", :castep_cell => ".cell", :castep_geom => ".cell",
                               :qe_in => ".in", :qe_out => ".in", :aims => "geometry.in", :xsf => ".xsf",
                               :pdb => ".pdb", :lammps => ".lmp", :car => ".car", :xsd => ".car", :res => ".res")

"""Save format matching a file (the .cell format for CASTEP .geom)."""
function save_extension(path::AbstractString)
    fmt = structure_format(path)
    haskey(SAVE_EXTENSION_OF, fmt) && return SAVE_EXTENSION_OF[fmt]
    e = lowercase(splitext(path)[2])
    return e in SAVE_EXTENSIONS ? e : ".vasp"
end

default_save_folder(v::Viewer) = v.path === nothing ? v.folder : dirname(v.path)

"""Whether Save may write straight back to `path` (see the notes at the top)."""
save_in_place(path) = path in WRITTEN || structure_format(path) in (:poscar, :xyz)

"""
    save_to!(viewer, path) -> Bool

Write the structure to `path` in the format given by its extension (or a
PNG image of the view for .png). A structure file becomes the viewer's file.
"""
function save_to!(v::Viewer, path::AbstractString)
    path = abspath(path)
    try
        if lowercase(splitext(path)[2]) == ".png"
            save_view(path, v)
            v.status[] = "Saved image $(basename(path)) in $(dirname(path))"
            return true
        end
        save_structure(path, v.crystal)
    catch err
        err isa InterruptException && rethrow()
        v.status[] = "Could not save $(basename(path)): " * sprint(showerror, err)
        return false
    end
    push!(WRITTEN, path)
    v.path = path
    v.saved[path] = v.crystal
    add_recent!(path)
    update_document!(v)
    v.status[] = "Saved $(basename(path)) in $(dirname(path))"
    return true
end

"""
    save!(viewer) -> Bool or the Save as tool

Save to the current file, or open Save as when there is none (or it should
not be overwritten).
"""
function save!(v::Viewer)
    v.path !== nothing && save_in_place(v.path) && return save_to!(v, v.path)
    return save_dialog(v)
end

"""Select every atom of the structure (Ctrl+A)."""
function select_all!(v::Viewer)
    v.selection = [(i, Shift(0, 0, 0)) for i in 1:natoms(v.crystal)]
    update_selection!(v)
end

# --- New structure -------------------------------------------------------------

"""1–3 lattice lengths: a (cubic), a c (hexagonal) or a b c."""
function parse_abc(s)
    w = tryparse.(Float64, split(strip(s), REAL_SEP; keepempty=false))
    (isempty(w) || length(w) > 3 || any(isnothing, w)) && return nothing
    x = Float64.(w)
    return length(x) == 1 ? (x[1], x[1], x[1]) : length(x) == 2 ? (x[1], x[1], x[2]) : Tuple(x)
end

"""
    new_structure_dialog(viewer) -> NamedTuple of the tool's controls

New structure tool: an empty cell (a b c, α β γ) to fill with Add, or a
prototype (diamond, zincblende, rocksalt, wurtzite, perovskite, ...) with
your elements and lattice constants. The new structure opens in a new tab
and is not saved until File → Save.
"""
function new_structure_dialog(v::Viewer)
    f = new_dialog(v, :new, "New structure")
    f === nothing && return nothing
    heading!(f, "New structure")
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], "Start from"; halign=:left, fontsize=PANEL_FONTSIZE)
    proto = Makie.Menu(f.g[r, 2]; options=prototype_names(), i_selected=1, fontsize=PANEL_FONTSIZE,
                       width=Makie.Relative(1))
    elements = field!(f, "Elements", "", s -> all(is_element, split(s)))
    abc = field!(f, "a b c (Å)", "", s -> parse_abc(s) !== nothing)
    angles = field!(f, "α β γ (°)", "90 90 90", is_real3)
    title = field!(f, "Title", "", _ -> true)
    hint = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], hint; halign=:left, fontsize=12, color=:gray40, justification=:left)
    create = button!(f, "Create")

    function show_prototype(k)
        name, nel, kind, els, lat, _ = PROTOTYPES[k]
        settext!(elements, els)
        settext!(abc, kind === :cubic ? string(lat[1]) : kind === :hexagonal ? "$(lat[1]) $(lat[3])" :
                                                                       join(lat, " "))
        hint[] = nel == 0 ? "An empty cell: add atoms with Add (panel) after Create." :
                 kind === :cubic ? "Cubic: type a. Elements in site order ($nel)." :
                 kind === :hexagonal ? "Hexagonal: type a c (γ = 120°). Elements in site order ($nel)." : ""
    end
    Makie.on(show_prototype, proto.i_selected)
    show_prototype(1)

    Makie.on(create.clicks) do _
        k = proto.i_selected[]
        name, nel, kind = PROTOTYPES[k]
        els = split(shown_text(elements))
        length(els) == nel && all(is_element, els) ||
            (v.status[] = nel == 0 ? "Enter no elements for an empty cell" :
                          "Enter $nel element" * (nel == 1 ? "" : "s") * " for $name, e.g. $(PROTOTYPES[k][4])"; return)
        lat = parse_abc(shown_text(abc))
        lat === nothing && (v.status[] = "Enter the lattice lengths a b c in Å"; return)
        ang = tryparse3(Float64, shown_text(angles))
        kind === :general && ang === nothing && (v.status[] = "Enter the angles α β γ in degrees"; return)
        create_structure!(v, name, els, lat, something(ang, [90.0, 90, 90]); title=shown_text(title))
    end
    finish_dialog!(v, :new, "New structure", f; labelwidth=0.3)
    return (; proto, elements, abc, angles, title, create)
end

"""A new, unsaved structure in a new tab (or in place of an untouched empty cell)."""
function create_structure!(v::Viewer, name, elements, abc, angles=(90, 90, 90); title::AbstractString="")
    c = try
        prototype_crystal(name, elements, abc, angles; title)
    catch err
        err isa InterruptException && rethrow()
        v.status[] = "Could not create $(lowercase(name)): " * sprint(showerror, err)
        return false
    end
    open_in_tab!(v, c)
    v.status[] = "New $(lowercase(name)): $(natoms(c)) atoms"
    return true
end

# --- Exit -------------------------------------------------------------------------

"""Unsaved changes: edits since the file was opened or saved (or any edit of an untitled structure)."""
is_modified(v::Viewer) = v.path === nothing ? !isempty(v.history) : get(v.saved, v.path, nothing) !== v.crystal

"""
    exit!(viewer; force=false) -> Bool

Close the AtomBond window (File → Exit, Ctrl+Q). With
unsaved changes the first Exit only warns in the status line; a second Exit
(with no edits in between) closes without saving. Returns whether it closed.
"""
function exit!(v::Viewer; force::Bool=false)
    unsaved = filter(k -> tab_modified(v, k), eachindex(v.tabs))
    state = [tab_crystal(v, k) for k in unsaved]
    if !force && !isempty(unsaved) && get(v.dialogs, :exit_armed, nothing) != state
        v.dialogs[:exit_armed] = state
        names = join((tab_label(v, k) for k in unsaved), ", ")
        v.status[] = "Unsaved changes in $names: choose Exit (Ctrl+Q) again to close without saving, or Save (Ctrl+S) first"
        return false
    end
    delete!(v.dialogs, :exit_armed)
    close_context_menu!(v)
    v.screen === nothing || (isopen(v.screen) && close(v.screen))
    v.status[] = "Closed"
    return true
end

# --- Menu bar -------------------------------------------------------------------

"""Items of menu `name` ("File", "Edit", "Build" or "View") for the current state."""
function menubar_items(v::Viewer, name::AbstractString)
    hassel = !isempty(v.selection)
    if name == "File"
        items = [MenuItem("New structure…", () -> new_structure_dialog(v); shortcut="Ctrl+N"),
                 MenuItem("Open…", () -> open_dialog(v); shortcut="Ctrl+O"),
                 MenuItem("Structure library…", () -> library_dialog(v); shortcut="Ctrl+L"),
                 MenuItem("Close tab", () -> close_tab!(v); shortcut="Ctrl+W")]
        recent = recent_files()
        if !isempty(recent)
            push!(items, SEPARATOR, menu_item("Recent files", nothing, false))
            for p in recent
                push!(items, MenuItem("  " * shorten(basename(p), 24), () -> open_structure!(v, p);
                                      shortcut=shorten(basename(dirname(p)), 14)))
            end
        end
        push!(items, SEPARATOR,
              MenuItem("Save", () -> save!(v); shortcut="Ctrl+S"),
              MenuItem("Save as…", () -> save_dialog(v); shortcut="Ctrl+Shift+S"),
              SEPARATOR,
              MenuItem("Export DFT inputs…", () -> export_dialog(v); shortcut="QE, ABINIT, SIESTA, DFTK, …"),
              MenuItem("Run DFTK…", () -> dftk_dialog(v); shortcut="DFTK.jl"),
              MenuItem("Save image…", () -> save_dialog(v; format=".png"); shortcut="PNG"),
              SEPARATOR,
              MenuItem("Exit", () -> exit!(v); shortcut="Ctrl+Q"))
        return items
    elseif name == "Edit"
        undo_label = isempty(v.history) ? "Undo" : "Undo " * shorten(v.history[end][2], 22)
        redo_label = isempty(v.future) ? "Redo" : "Redo " * shorten(v.future[end][2], 22)
        return [menu_item(undo_label, () -> undo!(v), !isempty(v.history); shortcut="Ctrl+Z"),
                menu_item(redo_label, () -> redo!(v), !isempty(v.future); shortcut="Ctrl+Y"),
                SEPARATOR,
                menu_item("Select all", () -> select_all!(v), natoms(v.crystal) > 0; shortcut="Ctrl+A"),
                menu_item("Clear selection", () -> clear_selection!(v), hassel; shortcut="Esc"),
                menu_item("Delete selected atoms", () -> delete_selected!(v), hassel; shortcut="Del"),
                SEPARATOR,
                MenuItem("Edit atoms…", () -> atoms_dialog(v))]
    elseif name == "Build"
        nout = length(atoms_outside_cell(v.crystal))
        return [MenuItem("New structure…", () -> new_structure_dialog(v); shortcut="Ctrl+N"),
                MenuItem("From library…", () -> library_dialog(v); shortcut="Ctrl+L"),
                MenuItem("Supercell…", () -> supercell_dialog(v)),
                MenuItem("Surface, superlattice…", () -> surface_dialog(v)),
                MenuItem("Nanotube…", () -> nanotube_dialog(v)),
                SEPARATOR,
                MenuItem("Cell and constraints…", () -> cell_dialog(v)),
                SEPARATOR,
                MenuItem("Shift atoms…", () -> shift_dialog(v)),
                MenuItem("Rotate atoms…", () -> rotate_dialog(v)),
                MenuItem("Oxygen bridge…", () -> bridge_dialog(v)),
                menu_item("Wrap atoms into cell", () -> wrap_atoms!(v), nout > 0;
                          shortcut=nout > 0 ? "$nout outside" : ""),
                SEPARATOR,
                MenuItem("Primitive cell", () -> build!(primitive_cell, v, "go to primitive cell")),
                MenuItem("Conventional cell", () -> build!(conventional_cell, v, "go to conventional cell")),
                MenuItem("Original structure", () -> build!(_ -> v.original, v, "go to original cell"))]
    elseif name == "View"
        all_lengths = v.settings.bond_lengths === :all
        return [MenuItem("Reset view", () -> reset_camera!(v));
                MenuItem("Turn view 90° counter-clockwise", () -> rotate_view!(v, 90); shortcut="↺");
                MenuItem("Turn view 90° clockwise", () -> rotate_view!(v, -90); shortcut="↻");
                SEPARATOR;
                MenuItem(all_lengths ? "Label selected atoms only" : "Label all atoms and bonds",
                         () -> set_bond_lengths!(v, all_lengths ? :selected : :all));
                MenuItem(v.settings.show_bond_lengths ? "Hide bond lengths" : "Show bond lengths",
                         () -> set_bond_length_labels!(v, !v.settings.show_bond_lengths); shortcut="Å");
                MenuItem(v.settings.show_bond_angles ? "Hide bond angles" : "Show bond angles",
                         () -> set_bond_angles!(v, !v.settings.show_bond_angles); shortcut="°");
                MenuItem(v.settings.show_dihedrals ? "Hide dihedral angles" : "Show dihedral angles",
                         () -> set_dihedrals!(v, !v.settings.show_dihedrals); shortcut="i–j–k–l");
                MenuItem(v.settings.show_charges ? "Hide atom charges" : "Show atom charges",
                         () -> set_charge_labels!(v, !v.settings.show_charges); shortcut="e");
                MenuItem(v.settings.show_magmoms ? "Hide magnetic moments" : "Show magnetic moments",
                         () -> set_magmom_labels!(v, !v.settings.show_magmoms); shortcut="μB");
                MenuItem(v.settings.show_coordinates ? "Hide atom coordinates" : "Show atom coordinates",
                         () -> set_coordinates!(v, !v.settings.show_coordinates); shortcut="x y z");
                MenuItem(coordinate_system(v.crystal, v.settings) === :fractional ? "Coordinates: use Cartesian (Å)" :
                         "Coordinates: use fractional",
                         () -> set_coordinates!(v, true; system=coordinate_system(v.crystal, v.settings) === :fractional ?
                                                               :cartesian : :fractional));
                SEPARATOR;
                [MenuItem("Along $label", () -> view_along!(v, which)) for (label, which) in VIEW_BUTTONS]]
    end
    throw(ArgumentError("no menu '$name'"))
end

const MENUBAR_TITLES = ["File", "Edit", "Build", "View"]

"""Open drop-down menu `k` of the menu bar under its title."""
function open_menubar_menu!(v::Viewer, k::Integer)
    b = v.dialogs[:menubar][k]
    bb = b.layoutobservables.computedbbox[]
    open_menu!(v, menubar_items(v, MENUBAR_TITLES[k]), (Makie.left(bb), Makie.bottom(bb));
               width=MENUBAR_MENU_WIDTH, owner=k)
end

"""Menu bar title at window pixel `xy`, or `nothing`."""
function menubar_title_at(v::Viewer, xy)
    haskey(v.dialogs, :menubar) || return nothing
    return findfirst(b -> xy in b.layoutobservables.computedbbox[], v.dialogs[:menubar])
end

function build_menubar!(pos, v::Viewer)
    g = Makie.GridLayout(pos; tellwidth=false, halign=:left, colgap=2)
    bar = map(enumerate(MENUBAR_TITLES)) do (k, title)
        b = Makie.Button(g[1, k]; label=title, fontsize=PANEL_FONTSIZE, padding=(10, 10, 4, 4),
                         buttoncolor=v.settings.background, buttoncolor_hover=Makie.RGBf(0.88, 0.91, 0.97),
                         buttoncolor_active=Makie.RGBf(0.80, 0.87, 0.98), strokewidth=0)
        Makie.on(_ -> open_menubar_menu!(v, k), b.clicks)
        b
    end
    v.dialogs[:menubar] = bar
    Makie.Label(g[1, length(bar)+1], v.document; fontsize=13, color=:gray35, halign=:right, tellwidth=false)
    Makie.colsize!(g, length(bar) + 1, Makie.Auto(false))
    return g
end
