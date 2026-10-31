# Supercells: new cell vectors are integer combinations of the old ones.

"""
    supercell(c, n1, n2, n3) -> Crystal
    supercell(c, P) -> Crystal

Repeat the cell n1×n2×n3 times, or more generally build the cell with
vectors (a', b', c') = (a, b, c) * P for an integer matrix P with det(P) > 0.
For example P = [1 -1 0; 1 1 0; 0 0 1] gives the √2×√2 rotated cell.
"""
supercell(c::Crystal, n1::Integer, n2::Integer, n3::Integer) = supercell(c, Diagonal([n1, n2, n3]))

function supercell(c::Crystal, P::AbstractMatrix{<:Integer})
    P = SMatrix{3,3,Int}(P)
    n = round(Int, det(P))
    n > 0 || throw(ArgumentError("supercell matrix must have a positive determinant (got $n)"))
    for k in 1:3
        # A non-periodic direction may not be mixed or repeated.
        e = SVector(ntuple(i -> Int(i == k), 3))
        c.pbc[k] || (P[:, k] == e && P[k, :] == e) ||
            throw(ArgumentError("cannot repeat along non-periodic direction $k"))
    end

    c = wrap(c)
    Pinv = inv(Float64.(P))
    # Old-cell translations that can land inside the new cell: the bounding
    # box of the new cell's corners in old fractional coordinates.
    corners = [P * SVector(i, j, k) for i in 0:1, j in 0:1, k in 0:1]
    lo = reduce((a, b) -> min.(a, b), corners)
    hi = reduce((a, b) -> max.(a, b), corners)
    ranges = [c.pbc[k] ? (lo[k]:hi[k]) : (0:0) for k in 1:3]

    species = Symbol[]
    frac = Vec3[]
    fixed = NTuple{3,Bool}[]
    charges, magmoms = Float64[], Float64[]
    eps = 1e-8
    for t in Iterators.product(ranges...), (s, f, x, q, m) in zip(c.species, c.frac, c.fixed, c.charges, c.magmoms)
        g = Pinv * (f + Vec3(t))
        # Inside the new cell along periodic directions (others are left as is).
        if all(k -> !c.pbc[k] || -eps <= g[k] < 1 - eps, 1:3)
            push!(species, s)
            push!(fixed, x)
            push!(charges, q)
            push!(magmoms, m)
            push!(frac, Vec3(ntuple(k -> c.pbc[k] ? max(g[k], 0.0) : g[k], 3)))
        end
    end
    length(frac) == n * natoms(c) ||
        error("supercell produced $(length(frac)) atoms, expected $(n * natoms(c))")
    # Atom constraints are repeated; cell constraints still mean a, b, c when
    # the new vectors are multiples of the old ones.
    return Crystal(Lattice(c.lattice.matrix * P), species, frac; pbc=c.pbc, title=c.title, fixed, charges, magmoms,
                   cell_constraints=isdiag(P) ? c.cell_constraints : FREE_CELL)
end
