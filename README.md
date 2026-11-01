# AtomBond.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://quantumworm.github.io/AtomBond.jl/dev/)
[![CI](https://github.com/QuantumWorm/AtomBond.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/QuantumWorm/AtomBond.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/QuantumWorm/AtomBond.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/QuantumWorm/AtomBond.jl)

AtomBond is a crystal-structure builder and interactive ball-and-stick viewer
written in Julia (GLMakie). It is meant for preparing ab-initio (DFT)
calculations: open or build a structure, edit it, set relaxation constraints,
and export inputs for open-source codes such as Quantum ESPRESSO, ABINIT and
SIESTA.

![AtomBond window showing a GaN wurtzite supercell](docs/src/assets/window.png)

- Reads and writes the common structure formats (CIF, XYZ, Quantum
  ESPRESSO, ABINIT, SIESTA, LAMMPS, XSF, PDB, …), with a built-in library of crystals (from the
  Crystallography Open Database) and molecules (experimental geometries
  from the literature).
- Builds supercells, slabs, superlattices and molecules in a box; edits
  atoms with undo; keeps several structures open in tabs.
- Labels bond lengths, angles, dihedrals, coordinates, charges and magnetic
  moments, and detects space groups (spglib).
- Writes DFT inputs with k-point meshes, fixed atoms and cell constraints.

## Installation

AtomBond needs Julia 1.10 or newer and OpenGL 3.3. It is not registered yet:

```julia
using Pkg
Pkg.add(url="https://github.com/QuantumWorm/AtomBond.jl")
```

## Usage

```julia
using AtomBond
visualize("structure.cif")           # a file in any supported format
visualize(load_library("GaN"))       # a structure from the built-in library
visualize()                          # an empty cell, to build a new structure
```

The same functions work without a window:

```julia
si = load_library("Silicon")
spacegroup(si)                                   # (227, "Fd-3m")
write_qe("si.in", supercell(si, 2, 2, 2); ecutwfc=40, pseudo="{El}.UPF")
```

See the [documentation](https://quantumworm.github.io/AtomBond.jl/dev/) for
the window, building and editing, file formats, DFT inputs and the API.

## Development

```bash
julia --project=. -e 'using Pkg; Pkg.test()'                  # tests
julia --project=docs -e 'using Pkg; Pkg.instantiate()'         # docs, once
julia --project=docs docs/make.jl                              # build docs/build
julia --project=. scripts/run_atombond.jl [file]               # open a window
```

`scripts/build_element_table.jl` and `scripts/build_library.jl` rebuild the
element table and the structure library from their published sources.

## License

MIT, see [LICENSE](LICENSE).
