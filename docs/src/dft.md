# DFT inputs and constraints

## Export

File → *Export DFT inputs…* writes:

- **VASP:** POSCAR and KPOINTS in a folder, plus `POTCAR.order` (the
  species order; POTCARs are licensed, so AtomBond does not write them).
  A structure with constraints also gets `INCAR.relax` (ISIF), and one with
  magnetic moments `INCAR.spin` (ISPIN, MAGMOM).
- **Quantum ESPRESSO:** a pw.x input, with `ecutwfc` and a pseudopotential
  file pattern such as `{El}.UPF` (`{El}` becomes each element symbol).
- **CASTEP:** a `.cell` file.
- **CIF** (P1) or **extended XYZ**.

The k-point mesh comes from a spacing in Å⁻¹ (the VASP KSPACING convention,
2π included); the tool shows the resulting mesh, and vacuum directions of
slabs and molecules get one k-point. In code: [`kpoint_grid`](@ref),
[`write_vasp_inputs`](@ref), [`write_qe`](@ref), [`write_castep`](@ref),
[`write_cif`](@ref), [`export_inputs`](@ref).

## Cell and constraints

Build → *Cell and constraints…* shows the cell (lengths, angles, volume,
density, vectors) and sets constraints for geometry optimization.

**Periodic box:** a molecule can be put into a periodic box with a given
vacuum (cubic or tight), so the cell is drawn and exported for plane-wave
codes; *Make it a molecule* goes back. In code: [`periodic_box`](@ref),
[`as_molecule`](@ref).

**Lattice constraints:** presets (all free, fixed cell, angles fixed, only c,
a = b and c, volume only) or numbers for a b c and α β γ in CASTEP's
convention (0 = fixed, equal numbers = vary together). Not every code can
express every combination:

| Code | How |
|------|-----|
| CASTEP | `CELL_CONSTRAINTS`, exact |
| Quantum ESPRESSO | `calculation = 'relax'`, or `'vc-relax'` with `cell_dofree` (`all`, `volume`, `xyz`, `z`) |
| VASP | `ISIF` 2 / 3 / 7 in `INCAR.relax`; other combinations need a patched VASP |

**Atom constraints:** fix the fractional x, y and/or z of the selected atoms
(also *Fix atoms* / *Free atoms* in the right-click menu). Fixed atoms are
marked × (+ when only some coordinates are fixed).

| Code | How |
|------|-----|
| VASP | POSCAR *Selective dynamics* |
| CASTEP | `IONIC_CONSTRAINTS` along the reciprocal directions |
| Quantum ESPRESSO | `if_pos` flags. pw.x applies them to Cartesian force components, which is exact unless a partly fixed axis is oblique; the file then says so |

POSCAR, CASTEP `.cell` and pw.x input read the constraints back. Constraints
stay through undo, shifts, rotations, wrapping and deletion, and supercells
repeat them. In code: [`fix_atoms`](@ref), [`free_atoms`](@ref),
[`with_constraints`](@ref), [`describe_cell_constraints`](@ref).

## DFTK and AtomsBase

[`to_atomsbase`](@ref) turns a structure into an AtomsBase system, which
DFTK and other AtomsBase-aware packages take directly;
[`from_atomsbase`](@ref) converts back. The conversion is tested; DFTK
itself is not a dependency and is not tested here.

```julia
using AtomBond, DFTK
system = to_atomsbase(load_library("Silicon"))
model  = model_DFT(system; functionals=PBE(), temperature=1e-3,
                   pseudopotentials=PseudoFamily("dojo.nc.sr.pbe.v0_4_1.standard.upf"))
```
