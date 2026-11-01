# Single-wall nanotubes as structures (periodic along the tube axis only).

include("nanotube_geometry.jl")

"""
Default C–C distance of a nanotube: that of graphite, a/√3 with a = 2.464 Å
(COD 9011577, Trucano & Chen, Nature 258, 136 (1975)).
"""
const GRAPHITE_CC = 2.464 / sqrt(3)

"""
    nanotube(n, m; bond=1.4226, periods=1, vacuum=10.0, species=(:C, :C)) -> Crystal

A single-wall (n, m) nanotube: the honeycomb sheet rolled along the chiral
vector n a₁ + m a₂. The tube axis is c, periodic, with `periods` translational
periods in the cell; a and b are square with `vacuum` Å between the walls of
neighbouring images and are not periodic (pbc = (false, false, true)).
`bond` is the nearest-neighbour distance before rolling (graphite's by
default); `species` puts one element on each sublattice, e.g. `(:B, :N)`
for a BN nanotube (give the B–N distance of hexagonal BN as `bond`). The
tube is an ideal rolled sheet: relax it before use.

    nanotube(6, 5)                 # chiral, semiconducting, 364 atoms
    nanotube(5, 5; periods=4)      # armchair, 4 periods along c
"""
function nanotube(n::Integer, m::Integer; bond::Real=GRAPHITE_CC, periods::Integer=1, vacuum::Real=10.0,
                  species=(:C, :C), title::AbstractString="")
    periods >= 1 || throw(ArgumentError("periods must be at least 1"))
    vacuum > 0 || throw(ArgumentError("vacuum must be positive"))
    length(species) == 2 || throw(ArgumentError("give two species, one per sublattice, e.g. (:C, :C) or (:B, :N)"))
    els = [known_element(s) for s in species]
    pos, sub, R, T = nanotube_geometry(n, m, bond)
    box = 2R + vacuum
    L = Lattice([box, 0, 0], [0, box, 0], [0, 0, periods * T])
    shift = [box / 2, box / 2, 0]
    cart = [Vec3(p + shift + [0, 0, k * T]) for k in 0:periods-1 for p in pos]
    spc = [els[s] for _ in 1:periods for s in sub]
    isempty(title) && (title = "($n,$m) " * (els[1] == els[2] == :C ? "carbon" : join(unique(els))) * " nanotube")
    return crystal_from_cartesian(L, spc, cart; pbc=(false, false, true), title)
end
