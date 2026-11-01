# DFT inputs and constraints

## Export

File → *Export DFT inputs…* writes:

- **Quantum ESPRESSO:** a pw.x input, with `ecutwfc` and a pseudopotential
  file pattern such as `{El}.UPF` (`{El}` becomes each element symbol).
- **ABINIT:** an `.abi` input (ABINIT 9 or newer) with `ngkpt`, `ecut` (the
  cutoff in Ha, half the Ry value entered) and `pseudos` from the pattern
  (`{El}.psp8` by default); magnetic moments give `nsppol 2` and `spinat`.
- **SIESTA:** an `.fdf` input with the k-grid, a DZP basis, PBE and
  `MeshCutoff 300 Ry` to adjust; pseudopotentials are read from
  `<Element>.psf` or `.psml`; magnetic moments give `DM.InitSpin`.
- **DFTK.jl:** a Julia script (`.jl`) for [DFTK](https://dftk.org) 0.8:
  the structure as an AtomsBase system, PBE with a little Gaussian smearing,
  `Ecut` (in Ha, half the Ry value entered), the k-grid and
  pseudopotentials from a PseudoPotentialData family (the PseudoDojo
  `dojo.nc.sr.pbe.v0_4_1.standard.upf` by default) or a `{El}` file
  pattern. It runs an SCF and prints energies and forces; with constraints
  it relaxes the structure with GeometryOptimization.jl instead. Magnetic
  moments go on the atoms (`magnetic_moment`). Run it with `julia name.jl`
  in an environment with DFTK, PseudoPotentialData, AtomsBase, Unitful and
  UnitfulAtomic (+ GeometryOptimization).
- **VASP:** POSCAR and KPOINTS in a folder, plus `POTCAR.order` (the
  species order; POTCARs are licensed, so AtomBond does not write them).
  A structure with constraints also gets `INCAR.relax` (ISIF), and one with
  magnetic moments `INCAR.spin` (ISPIN, MAGMOM).
- **CASTEP:** a `.cell` file.
- **CIF** (P1) or **extended XYZ**.

The written files have not been run through ABINIT or SIESTA by the tests
(they are read back by AtomBond's own readers); check them with your version
of the code. The DFTK scripts were run by hand with DFTK 0.8.1 (SCF, spin,
relaxation with a clamped atom); the test suite checks their contents
without installing DFTK.

The k-point mesh comes from a spacing in Å⁻¹ (2π included, n = ⌈|b|/Δk⌉);
the tool shows the resulting mesh, and vacuum directions of
slabs and molecules get one k-point. In code: [`kpoint_grid`](@ref),
[`write_qe`](@ref), [`write_abinit`](@ref), [`write_siesta`](@ref), [`write_dftk`](@ref),
[`write_vasp_inputs`](@ref), [`write_castep`](@ref), [`write_cif`](@ref),
[`export_inputs`](@ref).

## Cell and constraints

Build → *Cell and constraints…* shows the cell (lengths, angles, volume,
density, vectors) and sets constraints for geometry optimization.

**Periodic box:** a molecule can be put into a periodic box with a given
vacuum (cubic or tight), so the cell is drawn and exported for plane-wave
codes; *Make it a molecule* goes back. In code: [`periodic_box`](@ref),
[`as_molecule`](@ref).

**Lattice constraints:** presets (all free, fixed cell, angles fixed, only c,
a = b and c, volume only) or numbers for a b c and α β γ (0 = fixed, equal
numbers = vary together). Not every code can
express every combination:

| Code | How |
|------|-----|
| Quantum ESPRESSO | `calculation = 'relax'`, or `'vc-relax'` with `cell_dofree` (`all`, `volume`, `xyz`, `z`) |
| ABINIT | `ionmov 2` with `optcell` 0 (fixed), 1 (volume), 2 (all) or 6 (only c); other combinations are written as 2 with a note |
| SIESTA | `MD.TypeOfRun CG`, with `MD.VariableCell T` when the cell may change (the whole cell relaxes; other constraints are noted) |
| DFTK | `minimize_energy!` (GeometryOptimization) with `variablecell=true` when the cell may change (the whole cell relaxes; other constraints are noted) |
| VASP | `ISIF` 2 / 3 / 7 in `INCAR.relax`; other combinations need a patched VASP |
| CASTEP | `CELL_CONSTRAINTS`, exact |

**Atom constraints:** fix the fractional x, y and/or z of the selected atoms
(also *Fix atoms* / *Free atoms* in the right-click menu). Fixed atoms are
marked × (+ when only some coordinates are fixed).

| Code | How |
|------|-----|
| Quantum ESPRESSO | `if_pos` flags. pw.x applies them to Cartesian force components, which is exact unless a partly fixed axis is oblique; the file then says so |
| ABINIT | `iatfix` for fully fixed atoms; `iatfixx/y/z` for a fixed coordinate along a Cartesian axis (others are noted) |
| SIESTA | `Geometry.Constraints` for fully fixed atoms (partly fixed ones are noted) |
| DFTK | `clamp=true` on fully fixed atoms (partly fixed ones are noted) |
| VASP | POSCAR *Selective dynamics* |
| CASTEP | `IONIC_CONSTRAINTS` along the reciprocal directions |

POSCAR, CASTEP `.cell` and pw.x input read the constraints back. Constraints
stay through undo, shifts, rotations, wrapping and deletion, and supercells
repeat them. In code: [`fix_atoms`](@ref), [`free_atoms`](@ref),
[`with_constraints`](@ref), [`describe_cell_constraints`](@ref).

## DFTK and AtomsBase

Besides the exported script, DFTK can be used directly from Julia.
[`to_atomsbase`](@ref) turns a structure into an AtomsBase system, which
DFTK and other AtomsBase-aware packages take directly: magnetic moments
become the atoms' `magnetic_moment` and fully fixed atoms get `clamp=true`
(kept fixed by GeometryOptimization). [`from_atomsbase`](@ref) converts
back, keeping both. DFTK itself is not a dependency of AtomBond.

```julia
using AtomBond, DFTK, PseudoPotentialData
system = to_atomsbase(load_library("Silicon"))
model  = model_DFT(system; functionals=PBE(), temperature=1e-3,
                   pseudopotentials=PseudoFamily("dojo.nc.sr.pbe.v0_4_1.standard.upf"))
basis  = PlaneWaveBasis(model; Ecut=20, kgrid=collect(kpoint_grid(load_library("Silicon"))))
scfres = self_consistent_field(basis)

# relax, then look at the result in AtomBond
using GeometryOptimization
calc    = DFTKCalculator(; model_kwargs=(; functionals=PBE(), pseudopotentials=PseudoFamily("dojo.nc.sr.pbe.v0_4_1.standard.upf")),
                         basis_kwargs=(; Ecut=20, kgrid=[4, 4, 4]))
relaxed = minimize_energy!(system, calc).system
visualize(from_atomsbase(relaxed))
```
