# Using the window

## Layout

The window has a menu bar at the top, a row of tabs (one per open
structure), the structure's formula, space group and title, the 3D view, and
on the right either the view controls or one open tool. The line under the
view describes the selection; the line below it shows messages, in red when
something needs your attention.

Text boxes take effect when you click their button; pressing Enter first is
optional. Three numbers can be separated by spaces or commas, and supercell
repeats also by x (`2x2x1`).

## Mouse

- **Left-drag** rotates, the **scroll wheel** zooms, **right-drag** pans.
- **Click** an atom to select it, **Shift+click** to add more: 2 atoms show
  their distance, 3 the angle at the middle atom, 4 the dihedral, in click
  order. Click empty space or press **Esc** to clear. Atoms outside the home
  cell are marked with a prime (Si1′).
- **Right-click** an atom for its menu: *Edit atoms…*, *Change to X* for each
  other element in the structure, *Select all X*, *Delete*, *Fix / Free
  atoms*, *Shift / Rotate atoms…*, *Wrap atoms*, *Add oxygen bridge* (three
  atoms selected) and *Clear selection*. The clicked atom is selected first
  unless it is already part of the selection, so select several atoms and
  right-click one of them to act on all of them.

## Menus

| Menu | Items |
|------|-------|
| File | New structure… (Ctrl+N), Open… (Ctrl+O), Structure library… (Ctrl+L), Close tab (Ctrl+W), recent files, Save (Ctrl+S), Save as… (Ctrl+Shift+S), Export DFT inputs…, Save image…, Exit (Ctrl+Q) |
| Edit | Undo / Redo (named after the step), Select all (Ctrl+A), Clear selection, Delete selected atoms, Edit atoms… |
| Build | New structure…, From library…, Cell and constraints…, Supercell…, Surface / superlattice…, Shift atoms…, Rotate atoms…, Oxygen bridge…, Wrap atoms into cell, Primitive / Conventional cell, Original structure |
| View | Reset view, Turn view 90° counter-clockwise / clockwise, Label all atoms and bonds / selected atoms only, Show / Hide bond lengths, bond angles, dihedral angles, atom charges, magnetic moments, atom coordinates, Coordinates: fractional / Cartesian, view along a, b, c, a*, b*, c* |

The current file is shown at the right of the menu bar, with "modified"
after unsaved changes.

## View controls

The right-hand panel sets the display style (ball and stick, CPK, stick,
line), switches bonds, bond stubs (half-bonds to atoms outside the cell),
boundary atoms, the unit cell, atom labels and the axes, and sizes atoms and
bonds. *View along* aims the camera along a, b, c, a\*, b\*, c\* or any
[u v w]. The line of sight alone leaves the picture free to turn about it:
↺ 90° and ↻ 90° (also in the View menu) turn the view in the screen plane, and
*Turn* by any angle (positive is counter-clockwise), e.g. to put a chosen axis
horizontal; the axes in the corner turn with it. In code:
[`rotate_view!`](@ref)`(viewer, 90)`. *Orthographic* switches the projection. The *Cell* and *Edit selected
atoms* rows give quick access to supercells, primitive / conventional cells,
and replacing, deleting, moving or adding atoms.

## Tools

The menu items that end in "…" (New structure, Structure library, Cell and
constraints, Supercell, Surface, Shift atoms, Rotate atoms, Oxygen bridge,
Edit atoms, Export) open in the right-hand column, in place of the view
controls. Tools that change the structure say which tab they work on
("on GaN.cif") and follow you when you switch tabs. One tool is open at a
time; its *Close* button (or [`close_tool!`](@ref)) brings the view controls
back. Only the Open / Save file browser is a separate window.

## Tabs

Several structures can be open at once, one per tab. Open…, recent files,
New structure and the library each open a new tab (except in place of an
untouched empty cell, as `visualize()` starts); a file that is already open
is shown in its tab. Every tab keeps its own structure, *Original*, undo
history, selection, view and file, and edits, Undo and Save act on the tab
shown.

- Click a tab to show it; **Ctrl+Tab** / **Ctrl+Shift+Tab** go to the next /
  previous tab.
- **×** or **Ctrl+W** closes a tab. With unsaved changes (marked "*" on the
  tab) the first close only warns; close again to discard them.
- **+** starts a new structure.

In code: [`new_tab!`](@ref), [`switch_tab!`](@ref), [`close_tab!`](@ref),
[`ntabs`](@ref), [`active_tab`](@ref).

## Undo, saving and exit

Every change of the structure is one undoable step: **Ctrl+Z** undoes,
**Ctrl+Y** (or Ctrl+Shift+Z) redoes, per tab.

*Save* writes back to the current file when it is a POSCAR or XYZ file, or
one AtomBond saved. A CIF or CASTEP file opened from disk is not overwritten
(AtomBond writes only the cell and atoms and would drop the rest of such a
file), so Save opens *Save as* instead, as it does for an untitled structure.
*Save as* is a file browser with a file name and a format (see
[File formats](formats.md)), or a PNG image of the view. Replacing an
existing file asks for a second click. Opened and saved files go into the
recent-files list, kept in `~/.atombond/recent.txt` (or `$ATOMBOND_HOME`).

*Open* (Ctrl+O) lists folders and the files AtomBond can read; click a file
to see its formula and cell, double-click (or Enter) to open it. Arrow keys
and Page Up/Down move the selection, Backspace goes up a folder.

*Exit* (Ctrl+Q) closes the window. With unsaved changes in any tab the first
Exit only warns and names the tabs; Exit again closes without saving.

## Click feedback

Buttons flash after a click, a chosen menu item flashes before the menu
closes, and a click in the 3D view leaves a short ring at the cursor. A
click sound plays if `paplay` is installed. `click_feedback!(flash=false,
sound=false)` turns them off.
