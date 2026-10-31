# Build the documentation:  julia --project=docs docs/make.jl
# (the first time: julia --project=docs -e 'using Pkg; Pkg.instantiate()')

using Documenter
using AtomBond

DocMeta.setdocmeta!(AtomBond, :DocTestSetup, :(using AtomBond); recursive=true)

makedocs(;
    modules=[AtomBond],
    sitename="AtomBond.jl",
    authors="QuantumWorm",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://quantumworm.github.io/AtomBond.jl",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "Using the window" => "gui.md",
        "Building and editing" => "building.md",
        "Labels and measurements" => "labels.md",
        "File formats" => "formats.md",
        "DFT inputs and constraints" => "dft.md",
        "Library and element data" => "library.md",
        "Working from Julia code" => "scripting.md",
        "API reference" => [
            "Structures and files" => "api_structures.md",
            "Analysis and building" => "api_building.md",
            "Rendering and window" => "api_window.md",
        ],
    ],
    checkdocs=:exports,
)

get(ENV, "CI", "false") == "true" && deploydocs(; repo="github.com/QuantumWorm/AtomBond.jl", devbranch="main")
