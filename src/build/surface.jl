# Surfaces, slabs and superlattices.
#
# `surface_cell` re-expresses a bulk crystal in a cell whose a and b vectors
# lie in the (hkl) plane (same volume, same atoms). `slab` repeats that cell
# along c and adds vacuum; `stack_cells` / `superlattice` put such cells on top of
# each other.

"""
    surface_basis(hkl) -> P

Integer 3×3 matrix (columns = new cell vectors in units of the old ones) with
the first two columns spanning the (hkl) plane and the third satisfying
hkl ⋅ c3 = 1, so det(P) = ±1. Construction as in ASE's `surface`
(extended gcd); the in-plane pair is chosen to be as short as possible.
"""
function surface_basis(L::Lattice, hkl)
    h, k, l = hkl
    all(iszero, hkl) && throw(ArgumentError("(0 0 0) is not a plane"))
    g = gcd(h, k, l)
    h, k, l = h ÷ g, k ÷ g, l ÷ g
    if count(iszero, (h, k, l)) == 2
        # Low-index planes: keep the original axes, the normal axis last.
        c1, c2, c3 = h != 0 ? ((0, 1, 0), (0, 0, 1), (sign(h), 0, 0)) :
                     k != 0 ? ((0, 0, 1), (1, 0, 0), (0, sign(k), 0)) :
                              ((1, 0, 0), (0, 1, 0), (0, 0, sign(l)))
    else
        _, p, q = gcdx(k, l)               # p k + q l = gcd(k, l)
        a1, a2, a3 = cell_vectors(L)
        # Pick p, q (shifted along the solution line) so c1 and c2 are as
        # orthogonal as possible.
        k1 = dot(p * (k * a1 - h * a2) + q * (l * a1 - h * a3), l * a2 - k * a3)
        k2 = dot(l * (k * a1 - h * a2) - k * (l * a1 - h * a3), l * a2 - k * a3)
        if abs(k2) > 1e-10
            i = -round(Int, k1 / k2)
            p, q = p + i * l, q - i * k
        end
        _, a, b = gcdx(p * k + q * l, h)    # a (pk+ql) + b h = 1
        c1 = (p * k + q * l, -p * h, -q * h)
        c2 = (0, l, -k) .÷ gcd(l, k)
        c3 = (b, a * p, a * q)
    end
    return reduce_surface_basis(L, SVector{3,Int}(c1), SVector{3,Int}(c2), SVector{3,Int}(c3))
end

"""
Tidy a surface basis without changing the lattice it spans: Gauss-reduce the
in-plane pair (shortest vectors, 90° ≤ γ ≤ 120°), then remove as much of the
in-plane component of c3 as integer steps allow. The result is right-handed.
"""
function reduce_surface_basis(L::Lattice, u::SVector{3,Int}, v::SVector{3,Int}, w::SVector{3,Int})
    M = L.matrix
    len2(x) = sum(abs2, M * x)
    ip(x, y) = dot(M * x, M * y)
    while true
        len2(u) > len2(v) && ((u, v) = (v, u))
        m = round(Int, ip(u, v) / len2(u))
        m == 0 && break
        v -= m * u
    end
    ip(u, v) > 1e-10 && (v = -v)                 # obtuse (or right) angle between a and b
    # c3 -= integer combination of u, v closest to its in-plane part.
    G = [len2(u) ip(u, v); ip(u, v) len2(v)]
    x = round.(Int, G \ [ip(u, w), ip(v, w)])
    w -= x[1] * u + x[2] * v
    P = hcat(u, v, w)
    det(P) < 0 && (P = hcat(v, u, w))            # right-handed
    return SMatrix{3,3,Int}(P)
end

"""
Rotate the lattice so a ∥ x and b lies in the xy plane (c then has a
positive z component). Fractional coordinates are unchanged.
"""
function standard_orientation(c::Crystal)
    a, b, _ = cell_vectors(c.lattice)
    e1 = normalize(a)
    e3 = normalize(cross(a, b))
    e2 = cross(e3, e1)
    R = vcat(e1', e2', e3')
    with_positions(c, c.frac; lattice=Lattice(R * c.lattice.matrix))
end

"""Interplanar spacing of the planes spanned by a and b: thickness of one cell along the normal."""
layer_spacing(L::Lattice) = cell_volume(L) / norm(cross(cell_vectors(L)[1], cell_vectors(L)[2]))

"""
    surface_cell(c, hkl; shift=0) -> Crystal

Bulk cell equivalent to `c` with a, b in the (hkl) plane and c pointing out
of it, in standard orientation (surface normal along z). `shift` (fraction of
one layer spacing) moves the origin along c, which picks the termination of
slabs cut from it. The structure is still periodic bulk.
"""
function surface_cell(c::Crystal, hkl; shift::Real=0)
    all(c.pbc) || throw(ArgumentError("surface_cell needs a structure periodic in all directions"))
    s = supercell(c, surface_basis(c.lattice, hkl))
    frac = map(s.frac) do f
        z = mod(f[3] - shift, 1.0)
        Vec3(f[1], f[2], z > 1 - 1e-8 ? 0.0 : z)
    end
    title = isempty(c.title) ? "($(join(hkl, " "))) cell" : "$(c.title), ($(join(hkl, " "))) cell"
    return standard_orientation(Crystal(s.lattice, s.species, frac; title))
end

"""
    best_termination(c, hkl) -> shift

Try a cut between every pair of neighbouring atomic planes of the (hkl)
surface cell and return the `shift` that breaks the fewest bonds (ties: the
lowest shift). For Si(111) this picks the cut between bilayers (one dangling
bond per surface atom) rather than inside a bilayer (three).
"""
function best_termination(c::Crystal, hkl)
    sc = surface_cell(c, hkl)
    zs = sort(unique(round.([f[3] for f in sc.frac]; digits=4)))
    # Midpoints between consecutive planes, including the wrap-around gap.
    cuts = [(zs[k] + (k < length(zs) ? zs[k+1] : zs[1] + 1)) / 2 for k in eachindex(zs)]
    bulk = length(find_bonds(sc))
    best, best_broken = 0.0, typemax(Int)
    for s in sort(mod.(cuts, 1.0))
        sl = slab(c, hkl; layers=2, vacuum=10, shift=s)
        broken = 2bulk - length(find_bonds(sl))
        broken < best_broken && ((best, best_broken) = (s, broken))
    end
    return best
end

"""
    slab(c, hkl; layers=4, vacuum=15.0, shift=:auto) -> Crystal

Cut a slab of `layers` (hkl) repeat units from bulk crystal `c`. One layer is
one period of the surface cell along its normal (e.g. one cubic lattice
constant for Si (001) from the conventional cell). With `vacuum > 0` the c
vector is made normal to the surface and `vacuum` Å of empty space separate
the slab from its periodic image (the atoms are centred). `shift` chooses the
termination (see [`surface_cell`](@ref)); `:auto` uses the cut that breaks
the fewest bonds ([`best_termination`](@ref)).
"""
function slab(c::Crystal, hkl; layers::Integer=4, vacuum::Real=15.0, shift=:auto)
    layers >= 1 || throw(ArgumentError("layers must be ≥ 1"))
    vacuum >= 0 || throw(ArgumentError("vacuum must be ≥ 0"))
    shift === :auto && (shift = best_termination(c, hkl))
    s = supercell(surface_cell(c, hkl; shift), 1, 1, layers)
    s = Crystal(s.lattice, s.species, s.frac;
                title=replace(s.title, r"cell$" => "slab, $layers layers"))
    vacuum == 0 && return s
    return add_vacuum(s, vacuum)
end

"""
    add_vacuum(c, vacuum) -> Crystal

Replace c by a vector along the normal of the ab plane so that `vacuum` Å
of empty space lie between the top and bottom atoms of neighbouring images.
Atoms keep their Cartesian positions apart from a shift along the normal
that centres them; in-plane coordinates are wrapped into the cell.
"""
function add_vacuum(c::Crystal, vacuum::Real)
    a, b, _ = cell_vectors(c.lattice)
    n = normalize(cross(a, b))
    r = cartesian_positions(c)
    z = [dot(p, n) for p in r]
    zmin, zmax = extrema(z)
    Lz = zmax - zmin + vacuum
    Lz > 1e-6 || throw(ArgumentError("cell would have zero height"))
    L = Lattice(hcat(a, b, Lz * n))
    frac = map(r, z) do p, zi
        f = cart_to_frac(L, p)
        Vec3(mod(f[1], 1.0), mod(f[2], 1.0), (zi - zmin + vacuum / 2) / Lz)
    end
    with_positions(c, frac; lattice=L)
end

"""
    stack_cells(cells...; strain=false, tol=0.01) -> Crystal

Stack crystals along their c vectors (each block shifted by the c vectors of
the blocks below). The blocks must share the in-plane cell (a, b) within
relative tolerance `tol`; with `strain=true`, every block is strained in
plane to the first one (no Poisson relaxation along c).
"""
function stack_cells(cells::Crystal...; strain::Bool=false, tol::Real=0.01)
    isempty(cells) && throw(ArgumentError("nothing to stack"))
    ref = first(cells)
    a, b, _ = cell_vectors(ref.lattice)
    species = Symbol[]
    positions = Vec3[]
    offset = zero(Vec3)
    for (k, blk) in enumerate(cells)
        ak, bk, ck = cell_vectors(blk.lattice)
        mismatch = max(norm(ak - a) / norm(a), norm(bk - b) / norm(b))
        if mismatch > tol
            strain || throw(ArgumentError(@sprintf(
                "block %d does not match the in-plane cell of block 1 (%.1f%% mismatch); use strain=true", k, 100mismatch)))
        end
        Lk = Lattice(hcat(a, b, ck))       # strained to the first block's a, b
        for (s, f) in zip(blk.species, blk.frac)
            push!(species, s)
            push!(positions, frac_to_cart(Lk, f) + offset)
        end
        offset += ck
    end
    L = Lattice(hcat(a, b, offset))
    frac = [cart_to_frac(L, p) for p in positions]
    title = join((isempty(x.title) ? formula(x) : x.title for x in cells), " / ")
    return Crystal(L, species, frac; title)
end

"""
    superlattice(A, B, hkl; layers=(4, 4), strain=false) -> Crystal

A/B superlattice along (hkl): `layers[1]` (hkl) repeat units of A followed
by `layers[2]` of B, periodic in all directions.
"""
function superlattice(A::Crystal, B::Crystal, hkl; layers=(4, 4), strain::Bool=false)
    nA, nB = layers
    blockA = supercell(surface_cell(A, hkl), 1, 1, nA)
    blockB = supercell(surface_cell(B, hkl), 1, 1, nB)
    sl = stack_cells(blockA, blockB; strain)
    return Crystal(sl.lattice, sl.species, sl.frac;
                   title="$(formula(A)) / $(formula(B)) ($(join(hkl, " "))) superlattice, $nA+$nB layers")
end
