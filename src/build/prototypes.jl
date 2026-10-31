# Starting structures for File → New: an empty cell, or a common prototype
# filled with the user's elements and lattice parameters.

"""
(name, number of elements, lattice kind, default elements, default a b c,
fractional sites as (element number, position)). Lattice kind :cubic uses a
only, :hexagonal a and c (γ = 120°), :general a b c α β γ.
"""
const PROTOTYPES = [
    ("Empty cell", 0, :general, "", (5.0, 5.0, 5.0), Tuple{Int,NTuple{3,Float64}}[]),
    ("Simple cubic", 1, :cubic, "Po", (3.35, 3.35, 3.35), [(1, (0.0, 0.0, 0.0))]),
    ("Body-centred cubic", 1, :cubic, "Fe", (2.87, 2.87, 2.87), [(1, (0.0, 0.0, 0.0)), (1, (0.5, 0.5, 0.5))]),
    ("Face-centred cubic", 1, :cubic, "Cu", (3.615, 3.615, 3.615),
     [(1, (0.0, 0.0, 0.0)), (1, (0.0, 0.5, 0.5)), (1, (0.5, 0.0, 0.5)), (1, (0.5, 0.5, 0.0))]),
    ("Diamond", 1, :cubic, "Si", (5.431, 5.431, 5.431), :diamond),
    ("Zincblende", 2, :cubic, "Ga As", (5.653, 5.653, 5.653), :zincblende),
    ("Rocksalt", 2, :cubic, "Na Cl", (5.64, 5.64, 5.64), :rocksalt),
    ("Caesium chloride", 2, :cubic, "Cs Cl", (4.12, 4.12, 4.12), [(1, (0.0, 0.0, 0.0)), (2, (0.5, 0.5, 0.5))]),
    ("Fluorite", 2, :cubic, "Ca F", (5.463, 5.463, 5.463), :fluorite),
    ("Perovskite ABX₃", 3, :cubic, "Sr Ti O", (3.905, 3.905, 3.905),
     [(1, (0.0, 0.0, 0.0)), (2, (0.5, 0.5, 0.5)), (3, (0.5, 0.5, 0.0)), (3, (0.5, 0.0, 0.5)), (3, (0.0, 0.5, 0.5))]),
    ("Hexagonal close-packed", 1, :hexagonal, "Mg", (3.209, 3.209, 5.211),
     [(1, (1 / 3, 2 / 3, 0.25)), (1, (2 / 3, 1 / 3, 0.75))]),
    ("Wurtzite", 2, :hexagonal, "Ga N", (3.189, 3.189, 5.185),
     [(1, (1 / 3, 2 / 3, 0.0)), (1, (2 / 3, 1 / 3, 0.5)), (2, (1 / 3, 2 / 3, 0.375)), (2, (2 / 3, 1 / 3, 0.875))]),
]

const FCC_SITES = [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)]
shifted(sites, d) = [mod.(s .+ d, 1.0) for s in sites]

function prototype_sites(sites)
    sites === :diamond && return [(1, s) for s in [FCC_SITES; shifted(FCC_SITES, (0.25, 0.25, 0.25))]]
    sites === :zincblende && return [[(1, s) for s in FCC_SITES]; [(2, s) for s in shifted(FCC_SITES, (0.25, 0.25, 0.25))]]
    sites === :rocksalt && return [[(1, s) for s in FCC_SITES]; [(2, s) for s in shifted(FCC_SITES, (0.5, 0.0, 0.0))]]
    sites === :fluorite && return [[(1, s) for s in FCC_SITES]; [(2, s) for s in shifted(FCC_SITES, (0.25, 0.25, 0.25))];
                                   [(2, s) for s in shifted(FCC_SITES, (0.75, 0.75, 0.75))]]
    return sites
end

prototype_names() = [p[1] for p in PROTOTYPES]

"""
    prototype_crystal(name, elements, abc, angles=(90, 90, 90); title="") -> Crystal

A new structure: "Empty cell" (lattice `abc`, `angles`, no atoms) or one of
the prototypes in `PROTOTYPES` ("Diamond", "Zincblende", "Rocksalt",
"Wurtzite", ...) with `elements` (symbols, as many as the prototype needs) on
its sites. Cubic prototypes use `abc[1]` only; hexagonal ones `abc[1]` and `abc[3]`.
"""
function prototype_crystal(name::AbstractString, elements, abc, angles=(90, 90, 90); title::AbstractString="")
    k = findfirst(p -> p[1] == name, PROTOTYPES)
    k === nothing && throw(ArgumentError("unknown prototype '$name'"))
    _, nel, kind, _, _, sites = PROTOTYPES[k]
    els = [known_element(e) for e in elements]
    length(els) == nel || throw(ArgumentError("$name needs $nel element" * (nel == 1 ? "" : "s") *
                                              ", got $(length(els))"))
    all(>(0), abc) || throw(ArgumentError("lattice lengths must be positive"))
    L = kind === :cubic ? cubic(abc[1]) :
        kind === :hexagonal ? Lattice(abc[1], abc[1], abc[3], 90, 90, 120) :
        Lattice(abc..., angles...)
    s = prototype_sites(sites)
    isempty(title) && (title = nel == 0 ? "New structure" : "$(join(els)) $(lowercase(name))")
    return Crystal(L, Symbol[els[i] for (i, _) in s], [Vec3(p...) for (_, p) in s]; title)
end
