# DFT engines

File → *Run calculation…* runs the current structure with a DFT code, the
*engine*. [DFTK.jl](https://dftk.org) is the default: AtomBond installs it
into its own Julia environment the first time
([Running DFTK](dft.md#Running-DFTK-from-AtomBond)). AtomBond can also run
**Quantum ESPRESSO** (`pw.x`) and **CASTEP**, and you can add other codes
yourself (below).

AtomBond contains only the *adapter* for each code: the part that writes its
input, starts it and reads its output. The programs themselves, and their
pseudopotentials, are never part of AtomBond or its repository. They stay
on your machine, where AtomBond finds them as described below. This matters
for licensed codes such as CASTEP and VASP: everyone uses their own copy.

| Engine | Tasks | Results read back | Pseudopotentials |
|:-------|:------|:------------------|:-----------------|
| DFTK | SCF, band structure, relaxations, convergence check | everything (levels, DOS, bands, projections, steps) | PseudoDojo, downloaded automatically |
| Quantum ESPRESSO | SCF, relax, vc-relax (cell constraints as `cell_dofree`) | energy, forces, Fermi level, levels (DOS), every relaxation step with the stress | DFTK's PseudoDojo files (default), or a folder of UPF files |
| CASTEP | SinglePoint, GeometryOptimization (all cell constraints) | energy, forces, Fermi level, levels from `.bands` (DOS), steps from `.geom` with the stress | on the fly (default), or a folder |

All engines use PBE with a little Gaussian smearing (1 mHa), the cutoff
(Ry) and k-spacing of the form. Relaxations stop when the largest force is
below 1 meV/Å (DFTK, Quantum ESPRESSO) or 5 meV/Å (CASTEP), and, when the
cell relaxes, each stress component is below 0.05 GPa. The relaxation
window, the live progress and the DOS work the same way for all of them.

## Where AtomBond finds the programs

For each engine AtomBond uses the first of these:

1. **The settings file** `~/.config/atombond/engines.toml`. The *Save
   settings* button of the tool writes it; you can also edit it by hand.
   Each engine has a table with the full command (MPI included) and, if
   needed, the pseudopotential folder:

   ```toml
   [castep]
   command = "mpirun -np 8 castep.mpi"

   [qe]
   command = "mpirun -np 8 /opt/qe-7.4/bin/pw.x"
   pseudo_dir = "~/pseudo/sssp"
   ```

   On a shared machine or cluster, each user keeps their own file. Set
   `ATOMBOND_ENGINES` to use a file elsewhere, and `XDG_CONFIG_HOME` to
   move `~/.config`.
2. **An environment variable** `ATOMBOND_<ENGINE>_COMMAND`, e.g.
   `export ATOMBOND_CASTEP_COMMAND="mpirun -np 4 castep.mpi"`. This is
   useful in module scripts on clusters.
3. **The PATH**, under the usual program names: `pw.x`; `castep.serial`
   or `castep`. So a CASTEP copied into `~/bin` (which most Linux setups
   put on the PATH) is found without any setup.

The tool shows which one it found, e.g. *~/bin/castep.serial (from the
PATH)*. AtomBond adds the input file to the command (`pw.x -input Si.in`,
`castep Si`) and runs it in the job's folder, with the output in `Si.log`.

**Getting the open-source codes.** Quantum ESPRESSO is in most Linux
distributions and on conda-forge (`conda install -c conda-forge qe`, which
also has ABINIT and SIESTA). AtomBond does not install them.

**Pseudopotentials.** Quantum ESPRESSO reads UPF files. By default
(*PseudoDojo (as DFTK)* in the *Pseudopotentials* field) it uses the
PseudoDojo files that DFTK downloaded, so the two codes can be compared
directly and the tool's recommended cutoff applies. This needs DFTK
installed. To use other files, put a folder in the field, or set
`pseudo_dir` in the settings file or `$ESPRESSO_PSEUDO`; `~/espresso/pseudo`
is used when it exists. In a folder, each element's file is the one whose
name starts with the element symbol followed by `.`, `_` or `-`
(`Si.pbe-n-kjpaw_psl.1.0.0.UPF`, `si_pbe_v1.uspp.F.UPF`); the SSSP set
works too. CASTEP
generates its pseudopotentials on the fly unless the field names a folder
of `.usp`/`.recpot` files, which are then copied next to the input.

**Cutoff.** The recommended cutoff in the tool is for the PseudoDojo
pseudopotentials (DFTK, and Quantum ESPRESSO with the same files). With
other pseudopotentials, use the cutoff they need.

**Checked.** These runs used silicon from the library at 36 Ry with a
9×9×9 k-grid:

| Engine | Total energy | Gap on the k-grid | Relaxed a | Pressure before relaxing |
|:-------|:-------------|:------------------|:----------|:-------------------------|
| DFTK 0.8.1 | −230.2722 eV | 0.582 eV | 5.469 Å | +1.92 GPa |
| Quantum ESPRESSO 7.5 (same PseudoDojo files) | −230.2727 eV | 0.582 eV | 5.468 Å | +1.91 GPa |
| CASTEP 9.0 (on-the-fly pseudopotentials, 30 Ry) | (different pseudopotentials) | 0.588 eV | 5.469 Å | +1.95 GPa |

## Adding a DFT engine

An engine is a Julia type with a few methods: the input writer, the
program name and the output reader. To add one, for ABINIT, SIESTA, VASP
or an in-house code, write a file like the one below and `include` it
before starting AtomBond's window (or put it in
`~/.julia/config/startup.jl`). AtomBond already has input writers for
ABINIT, SIESTA and VASP ([`write_abinit`](@ref), [`write_siesta`](@ref),
`write_vasp_inputs`), and readers for many output formats
([`load_structure`](@ref)).

```julia
using AtomBond
import AtomBond: engine_id, engine_name, engine_executables, engine_command,
                 write_engine_inputs, read_engine_results, engine_tasks

struct ABINITEngine <: DFTEngine end

engine_id(::ABINITEngine) = :abinit                      # settings table [abinit], ATOMBOND_ABINIT_COMMAND
engine_name(::ABINITEngine) = "ABINIT"
engine_executables(::ABINITEngine) = ["abinit"]          # looked for on the PATH
engine_tasks(::ABINITEngine) = [:scf, :relax_atoms]      # what the Run tool offers

# Write the inputs into `dir`; return the main input and the file that exists when the run is done.
function write_engine_inputs(::ABINITEngine, dir, name, c::Crystal; relax, kspacing, ecutwfc, pseudo)
    input = joinpath(dir, name * ".abi")
    write_abinit(input, relax === :none ? c : fix_atoms(c, []); kspacing, ecutwfc, pseudo)
    return input, joinpath(dir, name * ".abo")
end

# How the program is started (in `dir`; its screen output goes to name.log).
engine_command(::ABINITEngine, command, dir, name, input) = `$command $(basename(input))`

# Read the finished run into a DFTResults (only the fields you have; the rest stay empty).
function read_engine_results(::ABINITEngine, job)
    text = read(job.results, String)
    m = match(r"etotal\s+([-\d.E+]+)", text)             # Ha
    return DFTResults(; calculation="scf", converged=true, energy=27.211386 * parse(Float64, m[1]),
                      structure=job.crystal, program="ABINIT")
end

register_engine!(ABINITEngine())
```

Optional methods:
- **[`engine_history`](@ref):** the steps of a running relaxation, for the
  live progress.
- **[`engine_finished`](@ref):** how to tell that the run really finished.
- **[`engine_pseudo_default`](@ref):** the default text of the
  pseudopotential field.
- **[`engine_install_hint`](@ref):** shown when the program is not found.

The fields of [`DFTResults`](@ref) are listed in its documentation. With
`eigenvalues`, `kweights` and `fermi` filled in, the DOS window works; with
`history` ([`RelaxHistory`](@ref)), the relaxation window works.

If the engine is useful to others, send it as a pull request to
`src/run/` with a test that uses a fake program (a shell script that
copies a prepared output file, as `test/runtests.jl` does for `pw.x` and
CASTEP). Never commit the program, its pseudopotentials, or output files
whose license forbids it.

## In code

[`run_dft`](@ref) runs any engine and returns the same kind of job as
[`run_dftk`](@ref):

```julia
using AtomBond
si = primitive_cell(load_library("Silicon"))
job = run_dft(si; engine=:castep, relax=:cell, ecutwfc=30)   # or engine=:qe, pseudo="~/pseudo/sssp"
wait(job)
r = job_results(job)           # DFTResults: energy, forces, levels, steps
plot_relaxation(r)
```

[`dft_engines`](@ref) lists the engines, [`find_engine`](@ref) says where a
program was found, and [`set_engine_setting!`](@ref) writes the settings
file.
