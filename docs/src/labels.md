# Labels and measurements

Selecting atoms measures them: 2 atoms give their distance, 3 the angle at
the middle atom, 4 the dihedral, shown on the structure and in the line
under the view. The View menu adds labels to the selected atoms, or to all
atoms and bonds with *Label all atoms and bonds* (the same item switches
back to *Label selected atoms only*).

| Label | View menu | Shown | With selected atoms |
|-------|-----------|-------|---------------------|
| Bond lengths | on by default (*Hide bond lengths*) | Å, 3 decimals, blue, at the bond's midpoint; stubs to atoms outside the cell too | 3 or more |
| Bond angles | *Show bond angles* | °, 1 decimal, green, inside the angle | 4 or more |
| Dihedral angles | *Show dihedral angles* | °, signed (IUPAC), purple, on a ring around the central bond | 5 or more |
| Coordinates | *Show atom coordinates* | (x, y, z), below-right of the atom | 1 or more |
| Charges | *Show atom charges* | e, red positive, blue negative, above-left | 1 or more |
| Magnetic moments | *Show magnetic moments* | μB, below-left | 1 or more |

Fewer selected atoms than listed show the plain measurement instead.
Torsions through a straight angle are left out, and −180° is shown as 180°.
Coordinates are fractional for crystals and Cartesian (Å) for molecules;
*Coordinates: use Cartesian / fractional* switches. A copy of an atom on a
cell face shows its own position (e.g. 1.000).

## Charges and magnetic moments

AtomBond does not compute or guess charges; they come from the file, and
atoms without one get no label.

- **Charges:** CASTEP `.castep` output (Mulliken, or Hirshfeld when there is
  no Mulliken table), LAMMPS data (charge and full styles), extended XYZ (a
  `charge`, `charges` or `initial_charges` column) and PDB formal charges.
  Extended XYZ and LAMMPS data save them.
- **Moments:** CASTEP `.cell` (`SPIN=`), CASTEP `.castep` output (Mulliken
  spins; one ħ/2 is 1 μB for g = 2) and extended XYZ (`magmom`, `magmoms`,
  `initial_magmoms`). They are saved to CASTEP `.cell`, extended XYZ, and by
  the VASP export as `INCAR.spin` with `ISPIN = 2` and `MAGMOM` in POSCAR
  order.

Supercells repeat charges and moments; changing an atom's element clears
them. Symmetry treats atoms with different moments or charges as different,
so an antiferromagnet keeps its magnetic cell: AFM bcc Cr is Pm-3m with 2
atoms, not Im-3m with 1 (`examples/structures/formats/Cr_afm.cell`).

## In code

[`set_bond_lengths!`](@ref)`(v, :all)`, [`set_bond_length_labels!`](@ref),
[`set_bond_angles!`](@ref), [`set_dihedrals!`](@ref),
[`set_coordinates!`](@ref), [`set_charge_labels!`](@ref),
[`set_magmom_labels!`](@ref), or the corresponding
[`ViewSettings`](@ref) fields, e.g.
`ViewSettings(bond_lengths=:all, show_bond_angles=true)`. For measurements
on coordinates: [`distance`](@ref), [`bond_angle`](@ref),
[`dihedral`](@ref), [`min_image_distance`](@ref).

## Bonds

Two atoms are bonded when their distance is at most 1.15 × the sum of their
covalent radii (at least the sum plus 0.15 Å, so H₂ is bonded). Metal–metal
bonds are left out when the structure also contains non-metals (otherwise
Ti–Ti in TiO₂ would show as bonds); pure metals and alloys keep them. The
limit can be set per element pair with `ViewSettings.bond_cutoffs`, or
`find_bonds(c; cutoffs)`.
