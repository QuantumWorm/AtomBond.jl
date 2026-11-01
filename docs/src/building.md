# Building and editing

All of these are undoable steps in the window, and functions that return a
new [`Crystal`](@ref) in code (see [Working from Julia code](scripting.md)).

## New structure

File → *New structure…* starts from an empty cell (a b c, α β γ; then add
atoms with *Add*) or a prototype filled with your elements and lattice
constants: simple, body-centred and face-centred cubic, diamond, zincblende,
rocksalt, CsCl, fluorite, perovskite ABX₃, hcp, wurtzite. Cubic prototypes
take a, hexagonal ones a and c. The new structure opens in a new tab and is
"Untitled" until saved. In code: [`prototype_crystal`](@ref).

## Editing atoms

The panel's *Edit selected atoms* row replaces the selected atoms by another
element, deletes them (also the Delete key), moves them by a Cartesian vector
in Å, or adds an atom at a fractional position. *Atoms…* opens the Edit
atoms tool: for one selected atom, change its element and set its exact
position (fractional or Cartesian); for several, change all their elements at
once. Its *Replace a species everywhere* section turns every atom of one
element into another (e.g. all As → P).

Changing elements keeps the bonds. After As → C in GaAs the Ga–C distance is
still 2.45 Å, longer than the usual Ga–C limit, so the new pair inherits the
old pair's limit and the status line says so (`ViewSettings.bond_cutoffs`).

In code: [`substitute`](@ref), [`delete_atoms`](@ref), [`move_atoms`](@ref),
[`add_atom`](@ref), [`set_position`](@ref), [`replace_species`](@ref).

## Supercells

Build → *Supercell…* makes repeats n1 n2 n3, or any cell whose vectors a′ b′
c′ are integer combinations of a, b, c (the columns of the supercell
matrix). The *Type* menu has presets: √2×√2 R45° (square ab), √3×√3 R30°
and rectangular 1×√3 (hexagonal ab). A preview shows the multiplicity, atom
count and new cell. *Primitive* and *Conventional* (panel or Build menu) use
spglib; *Original* goes back to the structure as opened. In code:
[`supercell`](@ref), [`primitive_cell`](@ref), [`conventional_cell`](@ref).

## Shift, rotate and wrap

**Shift atoms** moves all atoms, or the selected ones, by a vector in
fractional or Cartesian coordinates. *Selected atom → origin* fills in the
shift that puts the selected atom at (0, 0, 0). Atoms that end up outside
the cell are drawn where they are, not folded back in.

**Rotate atoms** turns all atoms, or the selected ones, by an angle about an
axis through a centre; the cell does not turn.

- Axis: a Cartesian direction, a lattice direction [u v w], the line from the
  first to the second selected atom (this also fixes the centre), or the
  viewing direction.
- Centre: the centroid of the rotated atoms, the first selected atom, or a
  point.
- Positive angles turn counterclockwise when the axis points toward you.
  Selected atoms turn as the copies you clicked, so a group that crosses a
  cell face turns as one piece.
- The preview names the closest contact afterwards and flags it below 0.6 ×
  the sum of covalent radii.

**Wrap atoms into cell** moves atoms outside the cell back in by whole
lattice vectors, so fractional coordinates end up in [0, 1); the periodic
structure does not change. Shift and Rotate have a switch to wrap afterwards.
`ViewSettings(wrap_atoms=true)` draws every atom inside the cell regardless.

In code: [`translate_atoms`](@ref), [`rotate_atoms`](@ref),
[`rotation_matrix`](@ref), [`closest_contact`](@ref), [`wrap`](@ref),
[`atoms_outside_cell`](@ref).

## Surfaces, slabs and superlattices

The Build surface tool (Build → *Surface, superlattice…*) cleaves a slab of
the current structure along (h k l) with a number of layers, vacuum and
termination (`auto` cuts where the fewest bonds break), or builds an A/B
superlattice of the current structure with a second structure file, with B
optionally strained in plane to A.

Superlattices are stacked from bulk blocks: interface spacings come from the
blocks and strained layers get no Poisson relaxation, so relax them.

In code: [`slab`](@ref), [`surface_cell`](@ref), [`add_vacuum`](@ref),
[`superlattice`](@ref), [`stack_cells`](@ref).

## Nanotubes

Build → *Nanotube…* makes a single-wall (n, m) nanotube for any indices:
n = m is armchair, m = 0 zigzag, anything else chiral, and a carbon tube is
metallic when n − m is a multiple of 3. The preview gives the diameter, the
period along the axis and the number of atoms before you create it; the tube
opens in a new tab. It is periodic along c only (pbc F F T), with the chosen
vacuum between neighbouring tubes in a and b, and as many periods along c as
asked. The bond length defaults to graphite's C–C distance (1.4226 Å, from
COD 9011577); two elements, e.g. `B N`, give a tube with one element per
sublattice (enter the B–N distance as the bond length). The tube is an ideal
rolled sheet, so relax it before use.

In code: [`nanotube`](@ref), e.g. `nanotube(6, 5)` or
`nanotube(5, 5; periods=4, vacuum=15)`.

## Bridging atoms

Build → *Oxygen bridge…* (or *Add oxygen bridge* in the right-click menu)
puts an atom into a bond. Select the two atoms of the bond, then a third atom
(Shift+click, in that order).

- The bridging atom (O by default; any element can be typed) goes to the
  middle of the bond plus a vector v of length *offset* (default 0.5 Å). v is
  perpendicular to the bond, in the plane of the three atoms, and points
  towards the third atom, so the new atom is inside their triangle.
- The bond between the two atoms is no longer drawn: they are now bonded
  through the bridge (the atoms themselves are not moved, so relax the
  structure to get realistic Si–O lengths and angles).
- *Other side of the bond (−v)* places it at −v instead; *Flip last bridge*
  moves a bridge you just added to the other side.
- The preview shows the two new bond lengths, the bridge angle, the closest
  other atom, and the largest offset that stays inside the triangle.

The bond is not stretched, so relax the result. With Si and the default
offset, Si–O is 1.28 Å and Si–O–Si 134°. In code: [`bridge_position`](@ref),
[`add_bridge`](@ref).
