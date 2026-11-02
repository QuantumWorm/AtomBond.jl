# API: running calculations

Running DFT calculations from AtomBond, with [DFTK.jl](https://dftk.org) or
another engine, and reading their results (see
[DFT inputs and constraints](dft.md) and [DFT engines](engines.md)). The plots of the
results are under [Rendering and window](api_window.md).

```@autodocs
Modules = [AtomBond]
Private = false
Pages = ["run/engines.jl", "run/dftk.jl", "run/qe.jl", "run/castep.jl"]
```
