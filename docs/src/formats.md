# File formats

[`load_structure`](@ref) and [`save_structure`](@ref) choose the format
from the file name.

| Code / program | Read | Write |
|---|---|---|
| Quantum ESPRESSO | pw.x input `.in`/`.pwi` (ibrav 0–12, -3, -12, 14; celldm or A…; `if_pos`); pw.x output `.out`/`.pwo` (last step) | pw.x input |
| ABINIT | `.abi`/`.abinit` (acell, rprim/angdeg, xred/xcart/xangst) | `.abi` input |
| SIESTA | `.fdf` | `.fdf` input |
| FHI-aims | `geometry.in` (`atom`, `atom_frac`) | `geometry.in` |
| VASP | POSCAR, CONTCAR, `*.vasp` (selective dynamics); XDATCAR (last frame); vasprun.xml (final structure) | POSCAR |
| CASTEP | `.cell` (with `SPIN=` and constraints); `.geom` (last frame); `.castep` output (final structure, charges, spins) | `.cell` |
| XCrySDen / VESTA | XSF (crystal, slab, polymer, molecule) | XSF |
| Gaussian and others | cube (atoms and cell; the grid is skipped) | – |
| Materials Studio | `.car`, `.xsd` | `.car` |
| SHELX / AIRSS | `.res`, `.ins` (with SYMM / LATT expansion) | `.res` |
| LAMMPS | data file (`.lmp`, `.data`, `data.*`; atomic / charge / full) | data file |
| General | CIF (symmetry expanded), XYZ / extended XYZ, PDB (CRYST1) | CIF (P1), extended XYZ, PDB |

- A `.in` file is read as pw.x, FHI-aims or ABINIT input depending on its
  contents.
- Files without a cell (XYZ, PDB without CRYST1, a molecule XSF, aims
  without lattice vectors) are put in a box with 5 Å of space and treated as
  molecules (not periodic). The box is drawn dashed.
- CIF: partial occupancies are drawn as full atoms, with a warning.
- Examples of each format are in `examples/structures/formats`. They are
  written by hand in each program's documented layout, not output of the
  programs themselves.

The separate input writers for DFT codes (with k-points and constraints) are
described in [DFT inputs and constraints](dft.md). [`to_atomsbase`](@ref) and
[`from_atomsbase`](@ref) convert to and from
[AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl) systems, which
packages such as DFTK accept.
