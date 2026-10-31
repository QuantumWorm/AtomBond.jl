# AtomBond

[![CI](https://github.com/QuantumWorm/AtomBond.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/QuantumWorm/AtomBond.jl/actions/workflows/CI.yml)

AtomBond is a crystal-structure builder and interactive ball-and-stick viewer
written in Julia (GLMakie). It is meant for preparing ab-initio (DFT)
calculations: open or build a structure, edit it, set relaxation constraints,
and export ready-to-run inputs for VASP, Quantum ESPRESSO or CASTEP. It works
much like the Materials Studio visualizer, but is free and open source.

![AtomBond window showing a GaN wurtzite supercell](docs/images/window.png)

## Features

- **Viewer:** ball and stick, CPK, stick and line styles; bonds across
  periodic boundaries; views along lattice and reciprocal axes or any [uvw];
  orthographic camera; PNG export.
- **Formats:** reads and writes CIF, VASP POSCAR, Quantum ESPRESSO, CASTEP,
  FHI-aims, XSF, XYZ, PDB, LAMMPS, Materials Studio `.car` and SHELX `.res`;
  also reads SIESTA, ABINIT, vasprun.xml, XDATCAR, cube and `.xsd`.
- **Building:** supercells (including general matrices), primitive and
  conventional cells, surfaces and slabs, superlattices, molecules in a
  periodic box, shifting, rotating and wrapping atoms, bridging atoms.
- **Editing:** select atoms with the mouse, then delete, substitute, add or
  move them, with undo and redo.
- **Analysis:** space group detection (spglib); labels for bond lengths,
  angles, dihedrals, coordinates, charges and magnetic moments.
- **DFT inputs:** fixed atoms and cell constraints, k-point meshes and spin,
  written to VASP, Quantum ESPRESSO and CASTEP input files.
- **Structure library:** 61 crystals and molecules from the Crystallography
  Open Database and NIST CCCBDB, each with its literature reference.
- **Element data:** masses, radii and colours from NIST, Cordero et al. and
  Jmol, in one table you can override.

![Methanol with bond lengths and angles](docs/images/methanol_labels.png)

## Installation

AtomBond needs Julia 1.10 or newer and a computer with OpenGL 3.3 (any
desktop; on Windows it also runs in WSL2). In Julia:

```julia
using Pkg
Pkg.add(url="https://github.com/QuantumWorm/AtomBond.jl")
```

Then open a structure:

```julia
using AtomBond
visualize("structure.cif")           # a file in any supported format
visualize(load_library("GaN"))       # or a structure from the built-in library
```

Left-drag rotates, the scroll wheel zooms and right-drag pans. From a script
(`julia script.jl`) keep the window open with `v = visualize(...); wait(v.screen)`.
The first start takes a minute or two while Julia compiles GLMakie.

## Running from a clone

```bash
# one-time: install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# open a structure in a window (left-drag rotate, scroll zoom, right-drag pan);
# the panel on the right switches styles, bonds/stubs, cell, views, builds
# supercells / primitive / conventional cells and saves PNGs
julia --project=. scripts/run_atombond.jl examples/structures/TiO2_rutile.cif

# run the tests
julia --project=. -e 'using Pkg; Pkg.test()'
```

Text boxes take effect when you click their button; pressing Enter first is
optional. Three numbers can be separated by spaces or commas, and supercell
repeats also by x (`2x2x1`). Messages that need your attention appear in red
at the bottom of the window.

In the window: **click** an atom to select it, **Shift+click** to add more
(2 atoms → distance, 3 → angle at the middle atom, 4 → dihedral, in click
order), click empty space or press **Esc** to clear. The result is shown on
the structure and in the line under it. Atoms outside the home cell are
marked with a prime (Si1′).

**Right-click** an atom to open a menu for it: it is selected first unless it
is already part of the selection, so select several atoms and right-click one
of them to act on all of them. The menu has *Edit atoms…* (the Edit atoms
window), *Change to X* for each other element in the structure, *Select all
X*, *Delete*, and *Clear selection*. Right-clicking empty space opens the same
menu when atoms are selected. A click outside the menu or **Esc** closes it,
and right-drag still pans.

**Editing** (panel section "Edit selected atoms"): type an element and press
*Replace* to substitute the selected atoms, *Delete* (or the Delete key) to
remove them, *Move* to shift them by a Cartesian vector in Å, or *Add* to put
a new atom of that element at a fractional position. Supercell / primitive /
conventional changes are edits too. *Atoms…* opens an Edit atoms window: for
one selected atom, change its element and set its exact position (fractional
or Cartesian Å; for an image atom, ′, the position of that copy); for several,
change all their elements at once. Its *Replace a species everywhere* section
turns every atom of one element into another (e.g. all As → P). **Ctrl+Z** undoes and **Ctrl+Y** (or
Ctrl+Shift+Z) redoes.

**Click feedback**: buttons flash blue after a click, the chosen menu item
flashes before the menu closes, and a click in the 3D view leaves a short
expanding ring at the cursor (deep orange on an atom, grey on empty space,
blue for a right-click). Effects are timed by drawn frames, so even a slow
first click shows them. A click sound plays if `paplay` is installed (WSLg
has an audio server; `sudo apt install pulseaudio-utils`).
`click_feedback!(flash=false, sound=false)` turns them off.

**Menu bar** (top of the window, like Materials Studio / VESTA):

| Menu | Items |
|------|-------|
| File | New structure… (Ctrl+N), Open… (Ctrl+O), Structure library… (Ctrl+L), recent files, Save (Ctrl+S), Save as… (Ctrl+Shift+S), Export DFT inputs…, Save image…, Exit (Ctrl+Q) |
| Edit | Undo / Redo (named after the step), Select all (Ctrl+A), Clear selection, Delete selected atoms, Edit atoms… |
| Build | New structure…, From library…, Cell and constraints…, Supercell…, Surface / superlattice…, Shift atoms…, Rotate atoms…, Oxygen bridge…, Wrap atoms into cell, Primitive / Conventional cell, Original structure |
| View | Reset view, Label all atoms and bonds / selected atoms only, Show / Hide bond lengths, bond angles, dihedral angles, atom charges, magnetic moments, atom coordinates, Coordinates: fractional / Cartesian, view along a, b, c, a*, b*, c* |

Moving the mouse to another title while a menu is open switches to it; a
click outside or Esc closes it. The current file is shown at the right of the
menu bar, with "modified" after unsaved changes (undoing back to the saved
state clears it). The panel's File row has New… / Open… / Save / Export….

**Supercell** (Build → *Supercell…*): repeats n1 n2 n3, or any cell whose
vectors a′ b′ c′ are integer combinations of a, b, c (columns of the
supercell matrix). The *Type* menu has presets: √2×√2 R45° (square ab),
√3×√3 R30° and rectangular 1×√3 (hexagonal ab); the repeats scale them. A
preview shows the multiplicity, atom count and new cell before you press
*Create supercell* (undoable). The panel's Supercell box still does plain
repeats.

**Shift atoms** (Build → *Shift atoms…*, or *Shift atoms…* in the right-click
menu): move all atoms, or only the selected ones, by a vector given in
fractional coordinates (of a, b, c) or Cartesian Å. Switching between the two
converts the typed vector. *Selected atom → origin* fills in the shift that puts
the selected atom at (0, 0, 0), e.g. to move the origin onto an atom. A
preview gives the shift in both coordinate systems and how many atoms will
end up outside the cell. Atoms outside the cell are drawn where they are (as
in Materials Studio), not folded back in.

**Rotate atoms** (Build → *Rotate atoms…*, or the right-click menu): turn all
atoms, or the selected ones, by an angle about an axis through a centre. The
cell itself does not turn.

- **Axis:**
  - a Cartesian direction x y z;
  - a lattice direction [u v w] (along u·a + v·b + w·c);
  - the line from the first selected atom you clicked to the second (e.g. turn
    a group about a bond; this also fixes the centre);
  - the viewing direction (turns the atoms in the screen plane).
- **Centre:** the centroid of the rotated atoms, the first selected atom, or a
  point in fractional or Cartesian coordinates.
- **Angle:** positive angles turn counterclockwise when the axis points toward
  you.
- **Groups across a cell face:** selected atoms turn as the copies you clicked,
  so a group that crosses a face turns as one piece.
- **Preview:** names the closest contact afterwards and flags it "too close!"
  below 0.6 × the sum of covalent radii. It also says how many atoms will be
  outside the cell, and warns when all atoms of a crystal turn: unless that is
  a symmetry operation, the crystal changes.
- **Wrap:** an optional switch wraps the atoms afterwards.
- **In code:** `rotate_atoms(c, angle, axis; center, atoms, images)`,
  `rotation_matrix`, `closest_contact`.

**Wrap atoms into cell** (Build menu, the Shift atoms window, or the
right-click menu for selected atoms outside the cell): moves atoms outside
the cell back in by whole lattice vectors, so fractional coordinates end up in
[0, 1). The periodic structure and its bonds don't change. "Wrap" is the usual
name for this (ASE `wrap()`, OVITO "Wrap at periodic boundaries", VMD
`pbc wrap`). The Build menu shows how many atoms are outside. Shift has a
*Wrap into the cell after shifting* switch to do both at once. In code:
`translate_atoms(c, v; fractional, atoms)`, `wrap(c; atoms)`,
`atoms_outside_cell(c)`; `ViewSettings(wrap_atoms=true)` draws every atom
inside the cell regardless.

**New structure** (File → *New structure…*): start from an empty cell (a b c,
α β γ; then add atoms with *Add*) or a prototype filled with your elements and
lattice constants: simple / body-centred / face-centred cubic, diamond,
zincblende, rocksalt, CsCl, fluorite, perovskite ABX₃, hcp, wurtzite. Cubic
prototypes take a; hexagonal ones a and c. Choosing a prototype fills in an
example (e.g. Ga As, 5.653 Å). The new structure replaces the current one as
an undoable step and is "Untitled" until saved.

**Saving**: *Save* writes back to the current file when it is a POSCAR or
XYZ file, or one AtomBond saved. A CIF or CASTEP file opened from disk is not
overwritten by Save (AtomBond writes only the cell and atoms and would drop
the rest of such a file), so Save opens *Save as* instead, as it does for an
untitled structure. *Save as* is the file browser with a file name and a
format: POSCAR (.vasp), CIF, XYZ, CASTEP .cell, Quantum ESPRESSO .in,
FHI-aims geometry.in, XSF, PDB, LAMMPS data, Materials Studio .car,
SHELX/AIRSS .res, or a PNG image of the view. Clicking a file copies its name; replacing an existing
file asks for a second click on Save. Opened and saved files go into the
recent-files list (kept in `~/.atombond/recent.txt`, or `$ATOMBOND_HOME`).

**Opening files** (*File → Open…* or **Ctrl+O**): a browser window lists
folders and the files AtomBond can read (see **File formats** below).
*All files* lists everything. Click a file to see its formula and
cell, double-click (or *Open* / Enter) to load it. Arrow keys and Page Up/Down
move the selection, Backspace goes up a folder, and the mouse wheel scrolls.
You can also type a folder path. Opening is undoable and the opened structure
becomes the *Original*. The browser starts in the folder of the current file.
WSL has no native file dialog, so this is drawn by AtomBond itself.

**File formats** (the format comes from the file name):

| Code / program | Read | Write |
|---|---|---|
| VASP | POSCAR, CONTCAR, `*.vasp`; XDATCAR (last frame); vasprun.xml (final structure) | POSCAR |
| CASTEP | `.cell` (with SPIN= moments); `.geom` (last frame); `.castep` output (final structure + Mulliken or Hirshfeld charges, Mulliken spins) | `.cell` (with SPIN=) |
| Quantum ESPRESSO | pw.x input `.in`/`.pwi` (ibrav 0–12, -3, -12, 14; celldm or A…); pw.x output `.out`/`.pwo` (last relax step) | pw.x input |
| FHI-aims | `geometry.in` (`atom`, `atom_frac`) | `geometry.in` |
| SIESTA | `.fdf` | – |
| ABINIT | `.abi`/`.abinit` (acell, rprim/angdeg, xred/xcart/xangst) | – |
| XCrySDen / VESTA | XSF (crystal, slab, polymer, molecule) | XSF |
| Gaussian / any code | cube (atoms and cell; the grid is skipped) | – |
| Materials Studio | `.car`, `.xsd` | `.car` |
| SHELX / AIRSS | `.res`, `.ins` (with SYMM / LATT expansion) | `.res` |
| LAMMPS | data file (`.lmp`, `.data`, `data.*`; atomic / charge / full) | data file (atomic) |
| General | CIF (symmetry expanded), XYZ / extended XYZ, PDB (CRYST1) | CIF (P1), extended XYZ, PDB |

A `.in` file is read as pw.x, FHI-aims or ABINIT input depending on its
contents. Files without a cell (XYZ, PDB without CRYST1, molecule XSF, aims
without lattice vectors) are put in a box with 5 Å of space and drawn as
non-periodic. Examples of each format are in `examples/structures/formats`
(hand-made in each program's documented layout).

**Surfaces and superlattices** (Cell → *Surface…* opens a second window):
cleave a slab of the current structure along (h k l) with a number of layers,
vacuum and termination (`auto` cuts where the fewest bonds break), or build an
A/B superlattice of the current structure with a second structure file, with
B optionally strained in plane to A. Results are undoable in the main window.

**Oxygen bridge** (Build → *Oxygen bridge…*, or *Add oxygen bridge* in the
right-click menu when three atoms are selected): select the two atoms of a
bond, then a third atom (Shift+click, in that order).
- **Placement:** the bridging atom (O by default; any element can be typed)
  goes to the middle of the bond plus a small vector v. v is perpendicular to
  the bond and lies in the plane of the three atoms, pointing towards the
  third atom, so the O is inside their triangle. Its length is the offset
  (default 0.5 Å).
- **Other side:** *Other side of the bond (−v)* places it at −v instead.
  *Flip last bridge* (also in the right-click menu) moves a bridge you just
  added to the other side.
- **Preview:** shows where the atom goes, the two bond distances, the bridge
  angle, the closest other atom, and the largest offset that stays inside the
  triangle.
- **Undo:** every step is undoable.
- **Relax afterwards:** the bond is not stretched, and its Si–Si line is still
  drawn through the O, so relax the result. With Si and the default 0.5 Å
  offset, Si–O is 1.28 Å and Si–O–Si 134°.

In code: `bridge_position(a, b, c; offset, flip)` and `add_bridge(crystal, a, b, c; element, offset, flip)`.

**Export** (File → *Export DFT inputs…*): VASP (POSCAR + KPOINTS), Quantum ESPRESSO
(`.in`, pw.x), CASTEP (`.cell`), CIF (P1) or extended XYZ. The k-point mesh
comes from a spacing in Å⁻¹ (VASP KSPACING convention, 2π included) and the
dialog shows the resulting mesh; vacuum directions of slabs get one k-point.
QE needs `ecutwfc` and a pseudopotential file pattern (`{El}.UPF`). POTCARs
are licensed, so VASP export writes the required species order to
`POTCAR.order` instead.

From the Julia REPL (`julia --project=.`):

```julia
using AtomBond
si = load_structure("examples/structures/Si_conv.vasp")
v = visualize(si)                               # interactive window with control panel
view_along!(v, :c); view_along!(v, [1, 1, 0])   # a, b, c, a*, b*, c* or [uvw]
v.settings.bond_stubs = false; redraw!(v)       # any setting, then redraw
save_view("si_view.png", v)                     # just the 3D view
snapshot("si.png", si; settings=ViewSettings(style=:cpk))   # PNG, no window
find_bonds(si)                                  # Bond(i, j, shift, length) list

rut = load_structure("examples/structures/TiO2_rutile.cif")   # symmetry-expanded
spacegroup(rut)                                 # (136, "P4_2/mnm")
primitive_cell(si); conventional_cell(si)       # via spglib
supercell(si, 2, 2, 1)                          # or supercell(si, [1 -1 0; 1 1 0; 0 0 1])
set_crystal!(v, supercell(si, 2, 2, 2))         # show it in an open window

select!(v, 5); select!(v, 1; add=true)          # same as click / Shift+click
v.info[]                                        # "Si5 – Si1   distance 2.3517 Å"
bond_angle(p, q, r); dihedral(p, q, r, s)       # on Cartesian points
min_image_distance(si, 1, 8)                    # shortest periodic distance

sige = substitute(si, [5], :Ge)                 # editing functions return new crystals
delete_atoms(si, [1]); move_atoms(si, [5], [0.1, 0, 0]); add_atom(si, :O, [0.5, 0.5, 0.5])
edit!(c -> substitute(c, [5], :Ge), v, "make SiGe")   # same, undoable in the window
undo!(v); redo!(v)
save_structure("SiGe.vasp", sige)

slab(si, (1, 1, 1); layers=4, vacuum=15)        # termination chosen automatically
surface_cell(si, (1, 1, 0))                     # bulk cell with a, b in the plane
ge = Crystal(cubic(5.658), fill(:Ge, 8), si.frac)
superlattice(si, ge, (0, 0, 1); layers=(4, 4), strain=true)   # Si/Ge (001) SL
stack_cells(cell1, cell2)                       # stack any cells sharing a, b

set_position(gaas, 5, (0.3, 0.25, 0.25))       # or (x, y, z) Å with cartesian=true
replace_species(gaas, :As, :P)                  # every As → P

kpoint_grid(si; kspacing=0.25)                  # (5, 5, 5)
write_qe("si.in", si; ecutwfc=40, pseudo="{El}.pbe-n-rrkjus_psl.1.0.0.UPF")
write_castep("si.cell", si); write_cif("si.cif", si)
write_vasp_inputs("si_vasp", si)                # POSCAR, KPOINTS, POTCAR.order
load_structure("relax.geom")                    # CASTEP: last frame of a .geom (or a .cell)
read_castep_geom_frames("relax.geom")           # every frame, titled with step and energy
export_inputs(si, :qe, "runs", "si")            # what the Export dialog calls
```

DFTK (or any AtomsBase-aware package) takes structures directly. The
conversion is tested; the DFTK calls below follow DFTK's documented API but
DFTK itself is not installed or tested here:

```julia
using AtomBond, DFTK
system = to_atomsbase(load_structure("examples/structures/Si_conv.vasp"))
model  = model_DFT(system; functionals=PBE(), temperature=1e-3,
                   pseudopotentials=PseudoFamily("dojo.nc.sr.pbe.v0_4_1.standard.upf"))
from_atomsbase(system)                          # and back to a Crystal
```

**Exit** (File → *Exit*, Ctrl+Q) closes the window and its dialog windows.
With unsaved changes the first Exit only warns in the status line; choosing
Exit again (without editing in between) closes without saving.

## Bond lengths, angles, dihedrals, coordinates, charges and moments

- **Selected atoms:** selecting 3 or more atoms labels the length (Å, 3
  decimals, blue) of every bond among them at the bond's midpoint. Two
  selected atoms show the measured distance as before.
- **All bonds:** View → *Label all atoms and bonds* labels every bond,
  including stubs to atoms outside the cell. The same item switches back
  to *Label selected atoms only*.
- **Bond angles:** View → *Show bond angles* adds the angles (°, 1 decimal,
  green) for the same bonds. Each label sits inside its angle, near the
  central atom. With selected atoms, angles need 4 or more selected, because
  3 show the measured angle. A 180° angle is placed beside the atom.
- **Dihedral angles:** View → *Show dihedral angles* adds the torsions
  i–j–k–l (°, signed, IUPAC convention, purple) about each bond j–k, for the
  same bonds. The labels sit on a ring around the bond, towards the two outer
  atoms, so the several dihedrals of one bond stay apart (anti pairs at right
  angles, nearer the bond). With selected atoms they need 5 or more, because
  4 show the measured dihedral. Torsions through a straight angle are left
  out, and −180° is shown as 180°.
- **Atom coordinates:** View → *Show atom coordinates* writes "(x, y, z)"
  (3 decimals, brown) below-right of each selected atom (one is enough), or
  of every atom with *Label all atoms and bonds*. Atom names, when on, go
  above-right. The coordinates are fractional for crystals and Cartesian Å
  for molecules; *Coordinates: use Cartesian / fractional* switches. A copy
  of an atom on a cell face shows its own position (e.g. 1.000).
- **Atom charges:** View → *Show atom charges* writes the charge in e
  ("+0.46", "−0.46"; red positive, blue negative) above-left of each atom.
  It uses the same atoms as the coordinates, and a single selected atom
  shows its charge in the info line too. Charges come from the file:
  - CASTEP `.castep` output: Mulliken charges, or Hirshfeld charges when
    there is no Mulliken table;
  - LAMMPS data (charge and full styles);
  - extended XYZ (a `charge`, `charges` or `initial_charges` column);
  - PDB formal charges (columns 79–80).

  Extended XYZ and LAMMPS (charge style) save them. Atoms without a known
  charge get no label; AtomBond does not guess oxidation states.
  Supercells repeat the charges, and changing an atom's element clears its
  charge.
- **Magnetic moments:** View → *Show magnetic moments* writes collinear
  moments in μB ("+2.15 μB" orange, "−2.15 μB" teal) below-left of the atoms
  (same atoms as the charges). A single selected atom shows "moment … μB"
  in the info line. Moments come from:
  - CASTEP `.cell` input (`SPIN=` on position lines);
  - CASTEP `.castep` output (Mulliken "Spin (hbar/2)" column; one hbar/2 is
    1 μB for g = 2);
  - extended XYZ (`magmom`, `magmoms`, `initial_magmoms` column).

  They are saved as:
  - CASTEP `.cell` `SPIN=`;
  - extended XYZ;
  - VASP export: an `INCAR.spin` with `ISPIN = 2` and `MAGMOM` in POSCAR
    order, e.g. `MAGMOM = 4*0 3*2 0`.

  Symmetry treats atoms with different moments (or charges) as different, so
  an antiferromagnet keeps its magnetic cell: AFM bcc Cr is Pm-3m with 2
  atoms, not Im-3m with 1. Example: `examples/structures/formats/Cr_afm.cell`.
- **Bond lengths off:** View → *Hide bond lengths* leaves only the other
  labels (e.g. charges).
- **In code:** `set_bond_lengths!(viewer, :all)`, `set_charge_labels!(viewer, true)`,
  `set_magmom_labels!(viewer, true)`, `set_bond_length_labels!(viewer, false)`,
  `set_bond_angles!(viewer, true)`, `set_dihedrals!(viewer, true)`,
  `set_coordinates!(viewer, true; system=:cartesian)`, or
  `ViewSettings(bond_lengths=:all, show_bond_angles=true, show_dihedrals=true)`.

## Cell and constraints

Build → *Cell and constraints…* has four sections:

- **Cell:** periodic or molecule, a b c α β γ, volume, density and the
  vectors.
- **Periodic box:** puts a molecule into a periodic box with the given vacuum
  (cubic or tight), so the cell is drawn and exported. *Make it a molecule*
  goes back.
- **Lattice constraints** for geometry optimization: presets (all free, fixed
  cell, angles fixed, only c, a = b and c, volume only) or numbers for a b c
  and α β γ in CASTEP's convention (0 = fixed, equal numbers = vary together).
  It shows how each code gets them:

  | Code | How |
  |------|-----|
  | CASTEP | `CELL_CONSTRAINTS`, exact |
  | QE | `calculation = 'relax'`, or `'vc-relax'` with `cell_dofree` (`all`, `volume`, `xyz`, `z`) |
  | VASP | `ISIF` 2 / 3 / 7, written to `INCAR.relax`; other combinations need a patched VASP |

- **Atom constraints:** fix the fractional x, y and/or z of the selected atoms
  (also *Fix atoms* / *Free atoms* in the right-click menu). Fixed atoms show
  × (+ when only some coordinates are fixed), and the info line says e.g.
  "fixed z".

How the exporters write atom constraints:

| Code | How |
|------|-----|
| VASP | POSCAR *Selective dynamics* (F = fixed; always fractional) |
| CASTEP | `IONIC_CONSTRAINTS` along the reciprocal directions |
| QE | `if_pos` flags. pw.x applies these to Cartesian force components, which is exact unless a partly fixed axis is oblique; the file then says so |

POSCAR, CASTEP `.cell` and pw.x input read the constraints back. Constraints
stay through undo, shifts, rotations, wrapping and deletion, and supercells
repeat them. In code: `fix_atoms(c, idx, (x, y, z))`, `free_atoms`,
`with_constraints(c; cell=(0,0,1,0,0,0))`, `periodic_box(c; vacuum, cubic)`.

## Structure library

File → *Structure library…* (Ctrl+L, also Build → *From library…*) lists
known crystals and molecules. It shows each one's formula, symmetry, cell and
source (database entry and paper), and *Open* loads it. Opening is undoable,
and the structure is Untitled until you save it, so Save never writes into the
library. In code: `structure_library()`, `load_library("Silicon")` (by name or
formula), `open_library!(viewer, "H2O")`.

| Group | Structures | Source |
|-------|------------|--------|
| Group IV | Si, diamond (COD 9011998, 9011997: Hom et al. 1975 precision lattice constants), Si0.5Ge0.5 (ordered model, a from Dismukes 1964, COD 1537802) | COD |
| III–V | GaAs, AlAs, InAs, GaP, AlP, InP, GaSb, AlSb, InSb (zincblende); GaN, AlN, InN (wurtzite) | 300 K lattice constants from Vurgaftman, Meyer, Ram-Mohan, *J. Appl. Phys.* 89, 5815 (2001); wurtzite u from Paszkowicz et al. (2003, 2004) Rietveld refinements |
| Oxides | α-quartz SiO2, corundum Al2O3, monoclinic HfO2 and ZrO2, rutile and anatase TiO2, ZnO, MgO, β-Ga2O3, SrTiO3 | COD (published refinements, listed per entry) |
| Diatomic molecules | H2, N2, O2, CO, NO, HF, HCl | NIST CCCBDB experimental geometries (all molecules) |
| Small molecules | H2O, NH3, CO2, H2S, N2O, SO2, O3, H2O2, HCN | |
| Organic molecules | CH4, C2H6, C2H4, C2H2, C6H6, CH3OH, C2H5OH, H2CO, HCOOH | |
| Process gases and precursors | SiH4, Si2H6, PH3, AsH3, B2H6, SiH2Cl2, SiCl4, SiF4, NF3, CF4, TiCl4 | |

- **Crystals:** CIFs from the [Crystallography Open Database](https://www.crystallography.net), an open
  collection of published experimental structures (it includes the American
  Mineralogist database). They are stored exactly as downloaded.
- **III–V semiconductors:** COD's III–V entries are mostly 1963 tables or
  older, and some are off (GaSb 6.118 Å, InAs 6.036 Å). So all twelve come
  from one review at 300 K, written as zincblende or wurtzite cells.
- **GaN:** a = 3.189 Å and c = 5.185 Å are bulk values. NIST's Structural
  Ceramics Database (Leszczynski 1994, bulk single crystal at 294 K) gives
  3.1880 / 5.18561 Å.
- **Molecules:** experimental geometries from NIST's
  [CCCBDB](https://cccbdb.nist.gov). Each has its point group and the papers CCCBDB cites. O2 and NO carry a note that their ground state is open-shell (use spin polarization). CCCBDB has no experimental geometry for GeH4, WF6, trimethylgallium or trimethylaluminium, so they are not included. Molecules open in a box with 5 Å around them, drawn as dashed grey lines (it is not periodic, but it is the cell written to input files).
  The window's *Molecules: in a periodic cubic box* switch (vacuum in Å)
  makes the box periodic instead, so the cell is shown and exported for
  CASTEP, VASP or QE, with 1 k-point along vacuum directions.
- **SiGe:** the real alloy is random (one site shared 50/50). The library
  holds the ordered model, with Si and Ge on the two diamond sublattices
  (F-43m), at the measured alloy lattice constant.

`data/library/` holds an index, `library.toml` (name, kind, category, formula,
file, description, source, url, reference, doi, note), with `crystals/` and
`molecules/` beside it. `scripts/build_library.jl` downloads the sources and
rebuilds the library. To add structures, add an entry there (COD id, or
CCCBDB CAS number) and run it again. Your own structures go into
`~/.atombond/library/library.toml`, in the same format with files next to it;
they are listed after the built-in ones.

## Element data

All element data is in one table, `data/elements.csv` (Z 1–103, one row per
element, one column per property, sources in its header). It is generated by
`scripts/build_element_table.jl` from the published sources:

| Column | Used for | Source |
|--------|----------|--------|
| `mass` | QE/LAMMPS masses, element from mass | NIST *Atomic Weights and Isotopic Compositions* (= IUPAC/CIAAW standard atomic weights; conventional value for interval elements; mass number of a long-lived isotope for radioactive ones) |
| `covalent_radius` | bonds (d ≤ 1.15 Σr), ball size (0.4 r, at least 0.2 Å) | Cordero et al., *Dalton Trans.* 2008, 2832 |
| `vdw_radius_bondi` | CPK size | Bondi, *J. Phys. Chem.* 68, 441 (1964), + Mantina et al. 2009 |
| `vdw_radius_alvarez` | CPK size where Bondi has none (e.g. transition metals, lanthanides) | Alvarez, *Dalton Trans.* 42, 8617 (2013) |
| `color_jmol` | atom colors | Jmol element colors (the CPK-style convention: H white, C grey, N blue, O red, S yellow, Cl green, P orange) |

No agency publishes radii for display. NIST and IUPAC cover masses; the
radii above are the standard literature sets that ASE, Jmol, OVITO and
VESTA-style tools use.

In code: `element(:Si)` gives an `Element` with `name`, `mass`,
`covalent_radius`, `vdw_radius` and `color`. `element_property(:Si,
:vdw_radius_alvarez)` reads any column, including ones added later.

**Your own values:** put an `elements.csv` with the same layout in
`~/.atombond/` (or `$ATOMBOND_HOME`). Give only the columns and rows to
change, e.g.

```
symbol,color,vdw_radius
O,FF2020,1.40
Si,1B3BFA,
```

`color` and `vdw_radius` take precedence over the built-in columns, and new
columns such as `electronegativity` become readable with `element_property`.
`use_vesta_colors!("…/VESTA/elements.ini")` takes the colors (and with
`radii=true` the van der Waals radii) from a VESTA installation. VESTA's
scheme is the other common one in papers, e.g. Si blue.

## Layout

| Folder | Contents |
|---|---|
| `src/core/` | `Element` data (read from `data/elements.csv`), `Lattice`, `Crystal` — no graphics |
| `data/elements.csv` | Element table: masses, radii, colors with their sources (built by `scripts/build_element_table.jl`) |
| `src/io/` | CIF reader (symmetry expansion) and P1 writer, POSCAR/CONTCAR and (extended) XYZ, CASTEP `.cell`/`.geom` readers, QE / CASTEP / VASP input writers, AtomsBase (DFTK) conversion, `load_structure` / `save_structure` |
| `src/analysis/` | periodic neighbor search (KD-tree), bond detection, space group / primitive / conventional cells (spglib), distances/angles/dihedrals |
| `src/build/` | supercells, editing (delete / substitute / add / move / set position / replace species), surfaces / slabs / superlattices, bridging atoms |
| `src/render/` | `ViewSettings` + styles, display model (what to draw), atoms, bonds, unit cell, camera presets, axes triad |
| `src/gui/` | `Viewer` window (`visualize`), control panel, mouse picking and selection, undo/redo + shortcuts, Build surface / Export / Open (file browser) / Edit atoms windows |
| `examples/structures/` | Si, GaAs, NaCl, GaN, water (POSCAR/XYZ); rutile TiO₂, ZnO, Si, GaN (CIF); GaAs (CASTEP .cell), a two-step Si relaxation (CASTEP .geom, hand-made) |
| `scripts/` | command-line launcher |
| `test/` | headless unit tests |

## Roadmap

0. ✅ Package skeleton, dependencies, GLMakie window under WSLg
1. ✅ Core model, element table, POSCAR/XYZ I/O, tests
2. ✅ Atoms as spheres, unit-cell box with colored a/b/c, legend, PNG export
3. ✅ Bonds across periodic boundaries (two-colored, half-bond stubs), boundary-atom images
4. ✅ Display styles (ball & stick, CPK, stick, line), control panel, view presets
   (a/b/c, a\*/b\*/c\*, [uvw]), orthographic, axes indicator, PNG export
5. ✅ CIF reading with symmetry expansion, space group detection, primitive/conventional cells, supercells
6. ✅ Click/Shift+click selection, atom labels, distance / angle / dihedral measurement
7. ✅ Editing (delete, substitute, add, move) with undo/redo; save POSCAR/XYZ
8. ✅ Surface cells, slabs (automatic termination), superlattices; Build surface dialog
9. ✅ Export: VASP, Quantum ESPRESSO, CASTEP, CIF, XYZ with k-point meshes; AtomsBase/DFTK conversion
10. ✅ File browser (Open…, Ctrl+O); CASTEP `.cell` and `.geom` reading
11. ✅ Edit atoms window (element / exact position of selected atoms, replace a species everywhere)
12. ✅ Right-click atom menu; menu bar (File / Edit / Build / View), New structure (prototypes), Save / Save as, recent files
13. ✅ Supercell dialog (repeats or general integer matrix, rotated-cell presets, live preview)
14. ✅ Click feedback (button / menu flash, ripple in the 3D view, optional click sound)
15. ✅ More formats: QE in/out, FHI-aims, SIESTA, ABINIT, XDATCAR, vasprun.xml, XSF, cube, PDB, LAMMPS, Materials Studio car/xsd, SHELX res
16. ✅ Shift atoms (fractional or Cartesian, all or selected) and Wrap atoms into cell
17. ✅ Rotate atoms (Cartesian / lattice / bond / view axis, centroid / atom / point centre, overlap check)
18. ✅ Element table from published sources (NIST masses, Cordero / Bondi / Alvarez radii, Jmol colors), user overrides, VESTA colors; File → Exit
19. ✅ Structure library: 25 crystals (group IV, III–V, oxides) and 36 molecules with sources
20. ✅ Cell and constraints: lattice / atom constraints for geometry optimization (VASP, CASTEP, QE), molecules in a periodic box
21. ✅ Labels for bond lengths, angles, dihedrals, coordinates, charges and magnetic moments
22. ✅ Oxygen bridge: an atom in a selected bond, inside or outside the triangle with a third atom

## Notes

- Under WSL2 the OpenGL renderer is Mesa *llvmpipe* (software). It works,
  but very large supercells will be slow until GPU passthrough is available.
- Bonds use covalent radii × 1.15. Metal–metal bonds are left out when the
  structure also contains non-metals (otherwise Ti–Ti in TiO₂, Hf–Hf in HfO₂
  would show as bonds); pure metals and alloys keep them.
- Changing elements keeps the bonds (as in Materials Studio). After As → C in
  GaAs the Ga–C distance is still 2.45 Å, longer than the Ga–C limit
  (1.15 × (1.22 + 0.76) = 2.28 Å). So the new pair inherits the old pair's
  limit (Ga–C ≤ 2.77 Å, as Ga–As), and the status line says so. The limit
  is per element pair (`ViewSettings.bond_cutoffs`, or `find_bonds(c;
  cutoffs)`), lasts for the session, and doesn't affect other pairs.
- CIF: partial occupancies are drawn as full atoms (with a warning).
- Superlattices are stacked from bulk blocks: interface spacings come from the
  blocks and strained layers get no Poisson relaxation, so relax them (DFT).
- Ambient occlusion (`ViewSettings(ssao=true)`) is off by default: on llvmpipe
  it only dims the whole picture. Try it on a machine with a real GPU.
