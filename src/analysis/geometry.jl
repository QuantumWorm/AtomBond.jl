# Distances, angles and dihedrals between Cartesian points (Å, degrees).

distance(p, q) = norm(q - p)

"""Angle p–q–r at the middle point q, in degrees."""
function bond_angle(p, q, r)
    u, v = p - q, r - q
    return acosd(clamp(dot(u, v) / (norm(u) * norm(v)), -1, 1))
end

"""
Dihedral (torsion) angle p–q–r–s in degrees, in (-180, 180]: the rotation
about q→r that takes the p–q–r plane onto the q–r–s plane (IUPAC sign).
"""
function dihedral(p, q, r, s)
    b1, b2, b3 = q - p, r - q, s - r
    n1, n2 = cross(b1, b2), cross(b2, b3)
    return atand(dot(cross(n1, n2), b2 / norm(b2)), dot(n1, n2))
end

"""
    measure(points) -> (kind, value) or nothing

Distance for 2 points, angle for 3 (at the middle one), dihedral for 4.
"""
function measure(points::AbstractVector)
    n = length(points)
    n == 2 && return (:distance, distance(points...))
    n == 3 && return (:angle, bond_angle(points...))
    n == 4 && return (:dihedral, dihedral(points...))
    return nothing
end

"""
    min_image_distance(c, i, j) -> Float64

Shortest distance between atom i and any periodic image of atom j.
"""
function min_image_distance(c::Crystal, i::Integer, j::Integer)
    d = c.frac[j] - c.frac[i]
    d = Vec3(ntuple(k -> c.pbc[k] ? d[k] - round(d[k]) : d[k], 3))
    # Rounding is exact for orthogonal cells; check neighbours for skewed ones.
    best = Inf
    for s in Iterators.product((c.pbc[k] ? (-1:1) : (0:0) for k in 1:3)...)
        best = min(best, norm(c.lattice.matrix * (d + Vec3(s))))
    end
    return best
end

"""Label used in the GUI, e.g. `"Si5"` (species + atom index)."""
atom_label(c::Crystal, i::Integer) = string(c.species[i], i)
