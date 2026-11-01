# Bond detection from covalent radii.

"""
A bond from atom `i` to the image of atom `j` displaced by `shift` lattice
vectors.
"""
struct Bond
    i::Int
    j::Int
    shift::Shift
    length::Float64
end

"""
    find_bonds(c; tolerance=1.15, min_length=0.4, metal_metal=:auto, cutoffs=Dict(), bridge_angle=80) -> Vector{Bond}

Atoms i and j are bonded when `min_length < d ≤ max(tolerance * (r_i + r_j), r_i + r_j + 0.3 Å)`,
with r the covalent radius. `cutoffs` can raise that limit for given element
pairs, e.g. `Dict((:C, :Ga) => 2.77)` (keys in alphabetical order, as `pairkey` makes them). Bonds
across periodic boundaries are included; each bond is listed once.

Covalent radii of metals are large enough that cation–cation pairs in oxides
and nitrides (Ti–Ti in rutile, Hf–Hf in HfO₂) would count as bonded. With
`metal_metal=:auto` such bonds are kept only when every species is a metal
(metals and alloys); pass `true`/`false` to force it.

A bond i–j is left out when an atom k bonded to both lies between them, i.e.
the angle i–k–j exceeds `bridge_angle` degrees: k bridges i and j (an O put
into a Si–Si bond, a bridging H). Three-membered rings keep all their bonds
(their angles are near 60°). `bridge_angle=nothing` keeps every bond.
"""
function find_bonds(c::Crystal; tolerance::Real=1.15, min_length::Real=0.4, metal_metal=:auto,
                    cutoffs::AbstractDict=Dict{Tuple{Symbol,Symbol},Float64}(), bridge_angle::Union{Nothing,Real}=80)
    natoms(c) == 0 && return Bond[]
    keep_mm = metal_metal === :auto ? all(ismetal, species_list(c)) : Bool(metal_metal)
    spc = species_list(c)
    limit = Dict(pairkey(a, b) => bond_cutoff(a, b; tolerance, cutoffs) for a in spc for b in spc)
    bonds = Bond[]
    for (i, j, s, d) in neighbor_pairs(c, maximum(values(limit)))
        si, sj = c.species[i], c.species[j]
        !keep_mm && ismetal(si) && ismetal(sj) && continue
        if min_length < d <= limit[pairkey(si, sj)]
            push!(bonds, Bond(i, j, s, d))
        end
    end
    return bridge_angle === nothing ? bonds : drop_bridged(c, bonds, bridge_angle)
end

"""
Bonds without those whose two atoms are bridged by a common neighbour k (the
angle i–k–j above `angle` degrees). Images are matched through the bond
vectors, so bridges across periodic boundaries count.
"""
function drop_bridged(c::Crystal, bonds::Vector{Bond}, angle::Real)
    bondvec(b) = frac_to_cart(c.lattice, c.frac[b.j] + Vec3(b.shift) - c.frac[b.i])
    nbrs = [Tuple{Int,Vec3}[] for _ in 1:natoms(c)]
    for b in bonds
        v = bondvec(b)
        push!(nbrs[b.i], (b.j, v)); push!(nbrs[b.j], (b.i, -v))
    end
    cosmax = cosd(angle)
    keep = trues(length(bonds))
    for (n, b) in enumerate(bonds)
        v = bondvec(b)                                  # i → j
        for (k, vik) in nbrs[b.i], (k2, vjk) in nbrs[b.j]
            k == k2 || continue
            norm(vik - (v + vjk)) < 1e-6 || continue    # the same image of k next to both
            a, bk = -vik, v - vik                       # k → i and k → j
            if dot(a, bk) / (norm(a) * norm(bk)) < cosmax
                keep[n] = false
                break
            end
        end
    end
    return bonds[keep]
end

"""Element pair in a fixed order, the key of bond `cutoffs`."""
pairkey(a::Symbol, b::Symbol) = a <= b ? (a, b) : (b, a)

"""
Longest distance at which elements `a` and `b` count as bonded: `tolerance`
times the sum of covalent radii, but at least the sum + 0.3 Å (which only
matters for light atoms: H–H in H2 is 0.741 Å, 1.15 × 0.62 = 0.713 Å; F–F in F2
is 1.414 Å, 1.15 × 1.14 = 1.311 Å).
"""
function bond_cutoff(a::Symbol, b::Symbol; tolerance::Real=1.15, cutoffs::AbstractDict=Dict())
    r = element(a).covalent_radius + element(b).covalent_radius
    return max(tolerance * r, r + MIN_BOND_SLACK, get(cutoffs, pairkey(a, b), 0.0))
end

const MIN_BOND_SLACK = 0.3      # Å

"""
    inherit_bond_cutoffs!(cutoffs, old, new; tolerance=1.15) -> pairs added

After a substitution (`new` has the atoms of `old` at the same places, some
with another element), let each element pair the substitution creates
between bonded atoms use at least the cutoff of the pair it replaced, so the
bonds stay (e.g. As → C in GaAs: Ga–C counts up to the Ga–As cutoff, 2.77 Å,
although Ga–C itself would stop at 2.28 Å). This is how Materials Studio
keeps bonds when an element is changed. Returns the pairs whose cutoff grew.
"""
function inherit_bond_cutoffs!(cutoffs::AbstractDict, old::Crystal, new::Crystal; tolerance::Real=1.15)
    natoms(old) == natoms(new) || throw(ArgumentError("not a substitution: atom counts differ"))
    added = Tuple{Symbol,Symbol}[]
    for b in find_bonds(old; tolerance, cutoffs)
        po = pairkey(old.species[b.i], old.species[b.j])
        pn = pairkey(new.species[b.i], new.species[b.j])
        po == pn && continue
        c_old = bond_cutoff(po...; tolerance, cutoffs)
        if c_old > bond_cutoff(pn...; tolerance, cutoffs)
            cutoffs[pn] = c_old
            push!(added, pn)
        end
    end
    return unique(added)
end

"""Number of bonds each atom takes part in."""
function coordination(c::Crystal, bonds::AbstractVector{Bond})
    z = zeros(Int, natoms(c))
    for b in bonds
        z[b.i] += 1
        z[b.j] += 1
    end
    return z
end
