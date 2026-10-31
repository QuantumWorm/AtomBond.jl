# Undoable edits and keyboard shortcuts.
#
# Every change of structure made through the window (edits, supercells,
# primitive/conventional cells) goes through `edit!`, which keeps the previous
# crystal on an undo stack.
#
#   Delete            delete the selected atoms
#   Ctrl+Z            undo
#   Ctrl+Y, Ctrl+Shift+Z   redo
#   Ctrl+A            select all atoms
#   Ctrl+N / Ctrl+O   new structure / open a file
#   Ctrl+S            save (Ctrl+Shift+S: save as)
#   Ctrl+L            structure library
#   Ctrl+Q            exit (close the window and its dialogs)
#
# Shortcuts are ignored while a text box in the panel has keyboard focus.

const MAX_HISTORY = 100

"""
    edit!(f, viewer, what; reframe=false, keep_selection=false)

Replace the crystal with `f(crystal)` as one undoable step. `what` describes
the step for the status line ("delete 2 atoms"). Errors from `f` are shown
in the status line and leave the structure unchanged. `reframe` re-aims the
camera (for changes of cell size); `keep_selection` keeps the selected atoms
(for edits that do not renumber them).
"""
function edit!(f, v::Viewer, what::AbstractString; reframe::Bool=false, keep_selection::Bool=false)
    c = try
        f(v.crystal)
    catch err
        err isa InterruptException && rethrow()
        v.status[] = "Could not $what: " * sprint(showerror, err)
        return false
    end
    # A pure substitution keeps the bonds of the substituted atoms drawn.
    kept = if natoms(c) == natoms(v.crystal) && c.frac == v.crystal.frac &&
              c.lattice.matrix == v.crystal.lattice.matrix && c.species != v.crystal.species
        vs = v.settings
        inherit_bond_cutoffs!(vs.bond_cutoffs, v.crystal, c; tolerance=vs.bond_tolerance)
    else
        Tuple{Symbol,Symbol}[]
    end
    push!(v.history, (v.crystal, what, reframe, v.path))
    length(v.history) > MAX_HISTORY && popfirst!(v.history)
    empty!(v.future)
    set_crystal!(v, c; reframe, keep_selection)
    v.status[] = "$(uppercasefirst(what)): $(natoms(c)) atoms   (Ctrl+Z undoes)" *
                 (isempty(kept) ? "" : "; bonds kept: " * join((@sprintf("%s–%s ≤ %.2f Å", p..., v.settings.bond_cutoffs[p])
                                                                 for p in kept), ", "))
    return true
end

function undo!(v::Viewer)
    isempty(v.history) && (v.status[] = "Nothing to undo"; return false)
    c, what, reframe, path = pop!(v.history)
    push!(v.future, (v.crystal, what, reframe, v.path))
    v.path = path                   # undoing an Open goes back to the previous file
    set_crystal!(v, c; reframe)
    v.status[] = "Undid: $what"
    return true
end

function redo!(v::Viewer)
    isempty(v.future) && (v.status[] = "Nothing to redo"; return false)
    c, what, reframe, path = pop!(v.future)
    push!(v.history, (v.crystal, what, reframe, v.path))
    v.path = path
    set_crystal!(v, c; reframe)
    v.status[] = "Redid: $what"
    return true
end

"""Crystal indices of the selected atoms (images count once), in click order."""
selected_indices(v::Viewer) = unique(first.(v.selection))

plural(n, word) = "$n $word" * (n == 1 ? "" : "s")

"""Run `f(indices)` on the selection, or say that nothing is selected."""
function with_selection(f, v::Viewer)
    idx = selected_indices(v)
    isempty(idx) ? (v.status[] = "Select atoms first (click / Shift+click)") : f(idx)
end

delete_selected!(v::Viewer) = with_selection(v) do idx
    edit!(c -> delete_atoms(c, idx), v, "delete $(plural(length(idx), "atom"))")
end

substitute_selected!(v::Viewer, el) = with_selection(v) do idx
    edit!(c -> substitute(c, idx, el), v, "change $(plural(length(idx), "atom")) to $el";
          keep_selection=true)
end

move_selected!(v::Viewer, shift) = with_selection(v) do idx
    edit!(c -> move_atoms(c, idx, shift), v,
          "move $(plural(length(idx), "atom")) by ($(join(shift, ", "))) Å"; keep_selection=true)
end

"""Three numbers as short text: "0.25 0 0"."""
vector_text(x) = join((@sprintf("%.6g", abs(t) < 1e-12 ? 0.0 : t) for t in x), " ")

outside_note(n) = n == 0 ? "" : "; $(plural(n, "atom")) " * (n == 1 ? "is" : "are") *
                                  " outside the cell (Build → Wrap atoms into cell)"

"""
    shift_atoms!(viewer, vector; fractional=false, atoms=nothing, wrap_after=false)

Undoable shift of all atoms (or the indices `atoms`) by `vector`, Cartesian
Å or fractional. `wrap_after` wraps the shifted atoms into the cell.
"""
function shift_atoms!(v::Viewer, vector; fractional::Bool=false, atoms=nothing, wrap_after::Bool=false)
    n = atoms === nothing ? natoms(v.crystal) : length(atoms)
    what = "shift $(plural(n, "atom")) by ($(vector_text(vector)))" * (fractional ? " (fractional)" : " Å") *
           (wrap_after ? " and wrap" : "")
    shifted(c) = (s = translate_atoms(c, vector; fractional, atoms); wrap_after ? wrap(s; atoms) : s)
    edit!(shifted, v, what; keep_selection=true) || return false
    v.status[] *= outside_note(length(atoms_outside_cell(v.crystal)))
    return true
end

"""
    wrap_atoms!(viewer; atoms=nothing)

Undoable wrap of the atoms outside the cell (or those of `atoms`) back into
it by whole lattice vectors.
"""
function wrap_atoms!(v::Viewer; atoms=nothing)
    out = atoms_outside_cell(v.crystal; atoms)
    isempty(out) && (v.status[] = "Nothing to wrap: the atoms are inside the cell"; return false)
    return edit!(c -> wrap(c; atoms), v, "wrap $(plural(length(out), "atom")) into the cell"; keep_selection=true)
end

"""
    rotate_atoms!(viewer, angle, axis; center, atoms, images, wrap_after=false, about)

Undoable rotation of all atoms (or `atoms`) by `angle` degrees about the
Cartesian `axis` through `center` (see [`rotate_atoms`](@ref)). `about`
names the axis in the step's description.
"""
function rotate_atoms!(v::Viewer, angle::Real, axis; center=nothing, atoms=nothing, images=nothing,
                       wrap_after::Bool=false, about::AbstractString="[$(vector_text(axis))]")
    n = atoms === nothing ? natoms(v.crystal) : length(atoms)
    what = "rotate $(plural(n, "atom")) by $(vector_text([angle]))° about $about" * (wrap_after ? " and wrap" : "")
    rotated(c) = (s = rotate_atoms(c, angle, axis; center, atoms, images); wrap_after ? wrap(s; atoms) : s)
    edit!(rotated, v, what; keep_selection=true) || return false
    v.status[] *= outside_note(length(atoms_outside_cell(v.crystal)))
    return true
end

# --- Constraints and periodicity ------------------------------------------------

axes_text(axes) = all(axes) ? "" : " (" * join(("x", "y", "z")[collect(Bool.(axes))], " ") * ")"

"""Fix fractional coordinates `axes` of the selected atoms (undoable)."""
fix_selected!(v::Viewer, axes=(true, true, true)) = with_selection(v) do idx
    any(axes) || (v.status[] = "Choose which coordinates to fix (x, y, z)"; return false)
    edit!(c -> fix_atoms(c, idx, axes), v, "fix $(plural(length(idx), "atom"))$(axes_text(axes))"; keep_selection=true)
end

"""Free the selected atoms (undoable)."""
free_selected!(v::Viewer) = with_selection(v) do idx
    edit!(c -> fix_atoms(c, idx, FREE_ATOM), v, "free $(plural(length(idx), "atom"))"; keep_selection=true)
end

"""Free every atom (undoable)."""
function free_all_atoms!(v::Viewer)
    n = count(any, v.crystal.fixed)
    n == 0 && (v.status[] = "Nothing to free: no atom is fixed"; return false)
    edit!(free_atoms, v, "free all $(plural(n, "fixed atom"))"; keep_selection=true)
end

"""Set the lattice constraints (a b c α β γ; 0 = fixed, equal = tied) (undoable)."""
function set_cell_constraints!(v::Viewer, cc)
    cc = try
        check_cell_constraints(Tuple(cc))
    catch err
        v.status[] = "Cannot use these cell constraints: " * sprint(showerror, err)
        return false
    end
    cc == v.crystal.cell_constraints && (v.status[] = "Nothing to change: the cell constraints are already " *
                                                      describe_cell_constraints(cc); return false)
    edit!(c -> with_constraints(c; cell=cc), v, "set cell constraints: $(describe_cell_constraints(cc))";
          keep_selection=true)
end

"""Put a molecule into a periodic box `vacuum` Å larger than it (undoable)."""
periodic_box!(v::Viewer; vacuum::Real=10.0, cubic::Bool=false) =
    edit!(c -> periodic_box(c; vacuum, cubic), v,
          "put the molecule in a periodic $(cubic ? "cubic " : "")box ($(vector_text([vacuum])) Å vacuum)";
          reframe=true, keep_selection=true)

"""Treat the structure as a molecule: not periodic, the cell only a container (undoable)."""
make_molecule!(v::Viewer) =
    !any(v.crystal.pbc) ? (v.status[] = "Nothing to change: the structure is already a molecule"; false) :
    edit!(as_molecule, v, "treat as a molecule (not periodic)"; reframe=true, keep_selection=true)

"""
    set_bond_lengths!(viewer, mode)

Bond-length labels: `:selected` (bonds among 3 or more selected atoms) or
`:all` (every bond).
"""
function set_bond_lengths!(v::Viewer, mode::Symbol)
    mode in (:selected, :all) || throw(ArgumentError("bond length labels are :selected or :all"))
    v.settings.bond_lengths = mode
    redraw!(v)
    v.status[] = mode === :all ? "Labelling all atoms and bonds (View menu to label only selected atoms)" :
                 "Labelling selected atoms (bonds: select 3 or more; angles: 4; dihedrals: 5)"
end

"""
    set_bond_angles!(viewer, on)

Show bond angles too, for the same bonds as the lengths (see
[`set_bond_lengths!`](@ref)).
"""
function set_bond_angles!(v::Viewer, on::Bool)
    v.settings.show_bond_angles = on
    redraw!(v)
    v.status[] = !on ? "Bond angles hidden" :
                 v.settings.bond_lengths === :all ? "Showing every bond angle" :
                 "Showing bond angles among selected atoms (select 4 or more; View → Label all atoms and bonds)"
end

"""
    set_dihedrals!(viewer, on)

Show dihedral angles too, for the same bonds as the lengths (see
[`set_bond_lengths!`](@ref)).
"""
function set_dihedrals!(v::Viewer, on::Bool)
    v.settings.show_dihedrals = on
    redraw!(v)
    v.status[] = !on ? "Dihedral angles hidden" :
                 v.settings.bond_lengths === :all ? "Showing every dihedral angle" :
                 "Showing dihedral angles among selected atoms (select 5 or more; View → Label all atoms and bonds)"
end

"""
    set_coordinates!(viewer, on; system=nothing)

Show atom coordinates beside the selected atoms (or all atoms, see
[`set_bond_lengths!`](@ref)). `system` = `:fractional`, `:cartesian` or
`:auto` changes the coordinate system.
"""
function set_coordinates!(v::Viewer, on::Bool; system::Union{Nothing,Symbol}=nothing)
    system === nothing || system in (:auto, :fractional, :cartesian) ||
        throw(ArgumentError("coordinates are :auto, :fractional or :cartesian"))
    v.settings.show_coordinates = on
    system === nothing || (v.settings.coordinates = system)
    redraw!(v)
    kind = coordinate_system(v.crystal, v.settings) === :fractional ? "fractional" : "Cartesian (Å)"
    v.status[] = !on ? "Atom coordinates hidden" :
                 v.settings.bond_lengths === :all ? "Showing $kind coordinates of every atom" :
                 "Showing $kind coordinates of the selected atoms (View → Label all atoms and bonds for all)"
end

"""
    set_charge_labels!(viewer, on)

Show atomic charges (e) beside the selected atoms (or all atoms, see
[`set_bond_lengths!`](@ref)). Charges come from the file: CASTEP .castep
(Mulliken), LAMMPS data (charge/full style), extended XYZ, PDB.
"""
function set_charge_labels!(v::Viewer, on::Bool)
    v.settings.show_charges = on
    redraw!(v)
    v.status[] = !on ? "Atom charges hidden" :
                 !has_charges(v.crystal) ? "Cannot show charges: this structure has none. They come from CASTEP .castep " *
                                           "output (Mulliken), LAMMPS data (charge/full), extended XYZ or PDB files" :
                 v.settings.bond_lengths === :all ? "Showing the charge of every atom (e)" :
                 "Showing the charges of the selected atoms (View → Label all atoms and bonds for all)"
end

"""
    set_bond_length_labels!(viewer, on)

Show or hide the bond-length labels (on by default), e.g. to see only
charges or angles.
"""
function set_bond_length_labels!(v::Viewer, on::Bool)
    v.settings.show_bond_lengths = on
    redraw!(v)
    v.status[] = on ? "Bond lengths shown" : "Bond lengths hidden"
end

"""
    set_magmom_labels!(viewer, on)

Show magnetic moments (μB) beside the selected atoms (or all atoms, see
[`set_bond_lengths!`](@ref)). Moments come from the file: CASTEP .castep
(Mulliken spins) and .cell (SPIN=), extended XYZ.
"""
function set_magmom_labels!(v::Viewer, on::Bool)
    v.settings.show_magmoms = on
    redraw!(v)
    v.status[] = !on ? "Magnetic moments hidden" :
                 !has_magmoms(v.crystal) ? "Cannot show magnetic moments: this structure has none. They come from " *
                                           "CASTEP .castep output (spin-polarized Mulliken), CASTEP .cell (SPIN=) " *
                                           "or extended XYZ (magmom column) files" :
                 v.settings.bond_lengths === :all ? "Showing the magnetic moment of every atom (μB)" :
                 "Showing the magnetic moments of the selected atoms (View → Label all atoms and bonds for all)"
end

# --- Bridging atoms -----------------------------------------------------------

"""
The three selected atoms in click order as (bond atom 1, bond atom 2,
plane/side atom), or a message saying what to select.
"""
function bridge_atoms(v::Viewer)
    sel = selected_atoms(v)
    length(sel) == 3 && return sel
    return "Select three atoms: the two bonded atoms, then a third for the plane and side (now $(length(sel)))"
end

"""
    add_bridge!(viewer; element=:O, offset=0.5, flip=false) -> Bool

Undoable: put an atom (O) into the bond between the first two selected atoms,
offset towards the third (inside their triangle) or away from it (`flip`).
"""
function add_bridge!(v::Viewer; element=:O, offset::Real=0.5, flip::Bool=false)
    sel = bridge_atoms(v)
    sel isa AbstractString && (v.status[] = sel; return false)
    a, b, third = (x.position for x in sel)
    el = known_element(element)
    names = join((display_label(v.display.crystal, x) for x in sel[1:2]), "–")
    ok = edit!(c -> add_bridge(c, a, b, third; element=el, offset, flip), v,
               "add $el bridge in $names ($(vector_text([offset])) Å $(flip ? "away from" : "towards") " *
               "$(display_label(v.display.crystal, sel[3])))"; keep_selection=true)
    ok && (v.dialogs[:bridge] = (; crystal=v.crystal, index=natoms(v.crystal), a, b, third, offset, flip, element=el))
    return ok
end

"""
    flip_bridge!(viewer) -> Bool

Undoable: move the bridging atom added last to the other side of its bond
(-v instead of v). Works while that structure is still the current one.
"""
function flip_bridge!(v::Viewer)
    last = get(v.dialogs, :bridge, nothing)
    (last === nothing || last.crystal !== v.crystal) &&
        (v.status[] = "Nothing to flip: add a bridge first (Flip works on the bridge just added)"; return false)
    flip = !last.flip
    pos = bridge_position(last.a, last.b, last.third; last.offset, flip)
    ok = edit!(c -> set_position(c, last.index, pos; cartesian=true), v,
               "flip the $(last.element) bridge to the other side of the bond"; keep_selection=true)
    ok && (v.dialogs[:bridge] = merge(last, (; crystal=v.crystal, flip)))
    return ok
end

function install_shortcuts!(v::Viewer)
    scene = v.lscene.scene
    Makie.on(Makie.events(scene).keyboardbutton; priority=10) do e
        e.action == Makie.Keyboard.press || return Makie.Consume(false)
        any(t -> t.focused[], v.textboxes) && return Makie.Consume(false)
        ctrl = Makie.ispressed(scene, Makie.Keyboard.left_control | Makie.Keyboard.right_control)
        shift = Makie.ispressed(scene, Makie.Keyboard.left_shift | Makie.Keyboard.right_shift)
        if e.key == Makie.Keyboard.delete
            delete_selected!(v)
        elseif ctrl && e.key == Makie.Keyboard.z
            shift ? redo!(v) : undo!(v)
        elseif ctrl && e.key == Makie.Keyboard.y
            redo!(v)
        elseif ctrl && e.key == Makie.Keyboard.o
            open_dialog(v)
        elseif ctrl && e.key == Makie.Keyboard.n
            new_structure_dialog(v)
        elseif ctrl && e.key == Makie.Keyboard.s
            shift ? save_dialog(v) : save!(v)
        elseif ctrl && e.key == Makie.Keyboard.a
            select_all!(v)
        elseif ctrl && e.key == Makie.Keyboard.l
            library_dialog(v)
        elseif ctrl && e.key == Makie.Keyboard.q
            exit!(v)
        else
            return Makie.Consume(false)
        end
        return Makie.Consume(true)
    end
end
