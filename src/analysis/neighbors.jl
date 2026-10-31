# Periodic neighbor search.
#
# All periodic images within reach of `cutoff` are put into a KD-tree; each
# atom of the home cell then queries the tree. Pairs come back once:
# (i, j, shift) with i < j, or i == j with a "positive" shift.

"""Distance between opposite faces of the cell along each lattice direction."""
function plane_spacings(L::Lattice)
    a, b, c = cell_vectors(L)
    V = cell_volume(L)
    return (V / norm(cross(b, c)), V / norm(cross(c, a)), V / norm(cross(a, b)))
end

"""How many cells to go out in each direction so every neighbor within `cutoff` is seen."""
function image_range(L::Lattice, cutoff::Real, pbc::NTuple{3,Bool})
    d = plane_spacings(L)
    ntuple(k -> pbc[k] ? ceil(Int, cutoff / d[k]) : 0, 3)
end

const Shift = SVector{3,Int}

ispositive(s::Shift) = s[1] > 0 || (s[1] == 0 && (s[2] > 0 || (s[2] == 0 && s[3] > 0)))

"""
    neighbor_pairs(c, cutoff) -> Vector{Tuple{Int,Int,Shift,Float64}}

All pairs (i, j, shift, distance) with |r_j + L*shift - r_i| ≤ cutoff, each
pair listed once.
"""
function neighbor_pairs(c::Crystal, cutoff::Real)
    n = natoms(c)
    pairs = Tuple{Int,Int,Shift,Float64}[]
    n == 0 && return pairs

    M = c.lattice.matrix
    home = cartesian_positions(c)
    nr = image_range(c.lattice, cutoff, c.pbc)
    shifts = vec([Shift(s1, s2, s3) for s1 in -nr[1]:nr[1], s2 in -nr[2]:nr[2], s3 in -nr[3]:nr[3]])

    points = Vector{Vec3}(undef, n * length(shifts))
    owner = Vector{Tuple{Int,Shift}}(undef, length(points))
    k = 0
    for s in shifts, j in 1:n
        k += 1
        points[k] = home[j] + M * s
        owner[k] = (j, s)
    end
    tree = KDTree(points)

    for i in 1:n
        for m in inrange(tree, home[i], cutoff)
            j, s = owner[m]
            (i < j || (i == j && ispositive(s))) || continue
            push!(pairs, (i, j, s, norm(points[m] - home[i])))
        end
    end
    return pairs
end
