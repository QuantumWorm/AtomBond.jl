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
  it relaxes the structure with GeometryOptimization.jl instead (in code,
  `relax=:atoms` or `:cell` asks for a relaxation of the atoms, or of the
  atoms and the cell, and `:none` for an SCF). Magnetic
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

## Running DFTK from AtomBond

File → *Run calculation…* runs the current structure in
[DFTK.jl](https://dftk.org) (PBE, 1 mHa Gaussian smearing), the default
engine; Quantum ESPRESSO, CASTEP or a code you add run from the same tool
(see [DFT engines](engines.md)). *Task*: an SCF
with or without the band structure, a relaxation of the atoms (*Relax
atoms*) or of the atoms and the cell (*Relax atoms and cell*, below), or a
convergence check (below); an SCF becomes a relaxation when the structure
has constraints (the same script as the DFTK export). The calculation runs in a
separate Julia process, so the window stays usable; the tool shows the last
lines of the log while it runs and, when it is done, a summary:

- whether the SCF (or relaxation) converged, and in how many iterations
  (steps);
- the total energy (eV and Ha) and the largest force (eV/Å);
- the Fermi level and the band gap: with the band structure, direct or
  indirect and where the band edges are (e.g. "indirect: VBM Γ, CBM
  between Γ and X"); without it, only the gap between the levels on the
  k-point grid, an estimate; a run with partly filled levels is reported
  as metallic.

**Band structure and dense DOS** (task *SCF + band structure*, the default): after the SCF
the script also computes, non-self-consistently from the SCF density, the
levels on a k-grid twice as dense, with more empty bands (for the DOS),
and, for a crystal without vacuum, the bands along the standard
high-symmetry path of its Brillouin zone (DFTK's `irrfbz_path`; for a
supercell, the path of the primitive cell). They go into `Name.bands.toml`;
a molecule or slab gets the dense DOS only. Before a band run the atoms
are put on exactly symmetric positions ([`symmetrize`](@ref), moves within
0.001 Å): DFTK and the k-path code use a much tighter symmetry tolerance
than AtomBond, and coordinates rounded in a CIF otherwise look like P1 to
them (α-quartz did; with the full symmetry it also ran twice as fast). After a relaxation the same
task first runs an SCF of the relaxed structure (with the symmetries it
has), so its levels, DOS and band structure are those of the result.

**Projected DOS.** Every SCF also projects its states on the pseudo-atomic
orbitals of the pseudopotentials (DFTK's orthonormalized
`atomic_orbital_projections`, the same idea as Quantum ESPRESSO's
`projwfc.x`), for the dense k-grid when the band step runs, into
`Name.pdos.toml`. The DOS and band-structure windows then draw the projected
DOS as coloured curves over the total DOS; *Projected DOS* chooses the
grouping: element and orbital ("Si 3s", "O 2p"), element, orbital type (s,
p, d), each atom (up to 12 atoms), or none. The curves add up to slightly
less than the total: the orbitals do not capture all of every state (the
"spilling", about 1 % for silicon). In code, [`pdos`](@ref) with `by` and
`atoms` (e.g. only the atoms around a defect).

**Relaxation steps.** A relaxation records every step in `Name.relax.toml`
(energy, the force on each atom, cell and positions). *Show relaxation*
opens a window with the energy and the largest force on an atom that may
move (clamped atoms keep their forces and are left out; DFTK's target,
1 meV/Å, dashed) against the step, and a *Step* slider: it shows that step
in a tab of the main window — *(DFTK relaxation, step k)* — with the forces
as red arrows (*Force arrows*; the largest force of the relaxation gets a
2 Å arrow, the same scale for every step). *Band structure / DOS of the
result* opens those plots when they were computed.

**Cell relaxation.** *Relax atoms and cell* also changes the cell, driven
by DFTK's stress, until the largest force is below 1 meV/Å and every stress
component below 0.05 GPa (Quantum ESPRESSO's default, 0.5 kbar); fixed
atoms stay clamped, but the whole cell relaxes (GeometryOptimization has no
partial cell constraints, so a strained-layer setup is not kept: use *Relax
atoms* for that). The atoms are first symmetrized, which keeps the stress
exactly symmetric. Each step records the stress too, and the window gets a
third panel with the volume change (%) and the pressure (GPa, dashed, zero
dotted); the summary gives the lattice parameters before and after, the
volume change and the remaining pressure. The cutoff matters more than for
energies: at a fixed cutoff the number of plane waves changes with the cell,
so the stress converges more slowly — use the *high* recommended cutoff or
check the lattice constant at two cutoffs. Silicon from the library
(a = 5.431 Å) at the recommended 36 Ry and a 9×9×9 k-grid: the pressure
starts at +1.92 GPa (the PBE cell wants to be larger) and one step later
it is −0.045 GPa, with a = 5.469 Å (+0.71 %, volume +2.15 %), in two SCF
runs. In code: [`pressure`](@ref), [`max_stress`](@ref),
[`cell_change`](@ref), [`variable_cell`](@ref).

While a relaxation runs, the tool's summary shows the latest step (energy
change, largest force) and *Show relaxation* already works: the window adds
the steps as they come and its slider moves on to the newest one — unless it
was moved back to an earlier step — and the relaxation tab of the main window
follows it while that tab is shown. In code: [`job_history`](@ref) reads the
steps of a running job; [`relax_figure`](@ref) and
[`plot_relaxation`](@ref) also take an `Observable{RelaxHistory}`; and
[`show_forces!`](@ref)
(arrows for any forces, e.g. those of an SCF).

*Open structure* opens the final structure (the relaxed one after a
relaxation) in a new tab. *Show DOS / bands* opens the band structure next
to the DOS (shared energy axis, band edges marked, gap in the subtitle;
*Save* writes `Name.bands.png`, `Name.bands.csv` and `Name.dos.csv`), or the
DOS alone when there is no band path, in its own window:

- every computed level is broadened into a Gaussian (the *Broadening*
  slider, 0.1 eV by default) and weighted with its k-point weight, in
  states per eV per cell; it integrates to the number of electrons below
  the Fermi level (dashed, at 0);
- filled states are shaded; spin down is drawn below the axis for a
  spin-polarized run;
- deep (semicore) levels, separated from the valence levels by more than
  5 eV, are hidden until *Show deep levels* is switched on;
- *Save* writes `Name.dos.png` and `Name.dos.csv` into the run folder.

Without the band step the DOS comes from the SCF levels only: few empty
bands, so it stops a little above the Fermi level. Every run gets its own folder,
`Folder/dftk_<Name>_<date-time>/`, with the script (`Name.jl`), the log
(`Name.log`) and the results (`Name.results.toml`: energies, forces,
eigenvalues and occupations per k-point, the final structure).

DFTK is large and is not a dependency of AtomBond: it is installed once into
its own Julia environment, by default `~/.julia/environments/atombond-dftk`.
The first *Run* does this by itself (DFTK, PseudoPotentialData, AtomsBase,
Unitful, UnitfulAtomic, GeometryOptimization; several minutes, mostly
compiling) and then starts the calculation, with the recommended cutoff
once the pseudopotentials are known; later runs start at once. *Install
DFTK* installs it without running anything. The PseudoDojo pseudopotentials are downloaded by
PseudoPotentialData on the first run.

**Cutoff.** The cutoff starts at the value the pseudopotentials recommend
for the elements of the structure, as soon as that is known (unless a
value was typed in). The cutoff line of the tool shows what the pseudopotential
family recommends for the elements of the structure — PseudoDojo gives a
low, normal and high cutoff per element; the structure needs the largest
of its elements — e.g. "Recommended 84 Ry (low 76, high 96)" for GaAs, and
warns when the cutoff is below the low one. *Use recommended* puts the
normal value in. The hints come from PseudoPotentialData in the DFTK
environment ([`dftk_cutoff_hints`](@ref), [`recommended_ecut`](@ref)); a
`{El}` file pattern has none.

**Convergence check.** Task *Convergence check* runs one SCF per cutoff (7
values from half to 1.5 times the recommended one, at the k-grid of the
form's k-spacing) and one per k-grid (k-spacings 0.6 to 0.15 1/Å, one per
distinct grid, at the form's cutoff; none for a molecule). Each total
energy is compared with the most accurate of its series; the converged
values are the smallest cutoff and largest k-spacing within 1 meV per atom
(for those and every more accurate step). The summary updates while it
runs; *Use converged values* puts them in the form and *Show plot* shows
|ΔE| per atom (log scale) and the gap against cutoff and k-spacing. The gap
on a k-grid jumps with the grid when the band edges lie between grid
points (silicon); the band structure gives the real gap. Silicon (PseudoDojo,
PBE): 32 Ry and a 9×9×9 grid (k-spacing 0.25) are within 1 meV/atom, in
14 SCF runs and about 2 minutes. Constraints are ignored in the check.

From Julia code:

```julia
using AtomBond
install_dftk() |> wait                     # once
job = run_dftk(primitive_cell(load_library("Diamond")); ecutwfc=40, bands=true)
job_status(job)                            # :running, then :done (or :failed, :stopped)
wait(job)
r = job_results(job)                       # DFTResults: energy, forces, levels, structure
r.energy, band_gap(r), r.structure
d = dos(r; width=0.1)                      # energies (eV), DOS per spin, Fermi level
write_dos("diamond.dos.csv", r)
plot_dos(r)                                # the DOS window
band_edges(r)                              # VBM, CBM, gap, direct?, where (from the band path)
write_bands("diamond.bands.csv", r)
plot_bands(r)                              # bands + DOS window
p = pdos(r; by=:element_orbital)           # projected DOS: "Si 3s" => …, "Si 3p" => …

c = primitive_cell(load_library("Silicon"))
recommended_ecut(c)                        # 36.0 Ry (PseudoDojo normal)
cj = run_dftk_convergence(c; kspacing=0.3) # cutoff and k-grid series
wait(cj)
cr = job_convergence(cj)
converged_settings(cr)                     # (ecut = 32.0, kspacing = 0.25)
plot_convergence(cr)
```

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
