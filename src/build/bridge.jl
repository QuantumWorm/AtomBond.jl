# Bridging atoms: put an O (or another atom) into the bond between two atoms.
#
# The bridge goes at the middle of bond a–b, moved by a small vector v that is
# perpendicular to the bond and lies in the plane of a, b and a third atom c,
# pointing towards c (so the bridge is inside triangle a b c); `flip` uses -v,
# the other side of the bond.

"""
    bridge_offset_direction(a, b, c) -> unit vector

Unit vector perpendicular to bond a–b, in the plane of a, b and c, pointing
towards c. Throws when c is on the line through a and b (no plane).
"""
function bridge_offset_direction(a, b, c)
    u = normalize(Vec3(b - a))
    w = Vec3(c - (a + b) / 2)
    w -= dot(w, u) * u
    norm(w) < 1e-6 * norm(b - a) &&
        throw(ArgumentError("the third atom is on the line of the bond, so it does not define a plane"))
    return normalize(w)
end

"""
    bridge_position(a, b, c; offset=0.5, flip=false) -> Vec3

Cartesian position of an atom bridging bond a–b: the bond's midpoint plus
`offset` Å along [`bridge_offset_direction`](@ref) (towards c, inside the
triangle a b c), or minus it with `flip=true` (the other side of the bond).
"""
function bridge_position(a, b, c; offset::Real=0.5, flip::Bool=false)
    offset >= 0 || throw(ArgumentError("the offset must be ≥ 0"))
    v = offset * bridge_offset_direction(a, b, c)
    return (a + b) / 2 + (flip ? -v : v)
end

"""
Largest offset that keeps the bridge inside triangle a b c (towards c).
In the plane, with the bond from -L/2 to +L/2 along x and c at (s, h), the
line x = 0 leaves the triangle through edge a–c (s ≥ 0) or b–c (s < 0), at
height h·(L/2)/(L/2 + |s|).
"""
function bridge_max_offset(a, b, c)
    m, d, u = (a + b) / 2, bridge_offset_direction(a, b, c), normalize(Vec3(b - a))
    h, s, half = dot(Vec3(c - m), d), dot(Vec3(c - m), u), norm(b - a) / 2
    return h * half / (half + abs(s))
end

"""
    add_bridge(c, a, b, third; element=:O, offset=0.5, flip=false, min_distance=0.5) -> Crystal

Add a bridging atom (oxygen by default) between Cartesian positions `a` and
`b` (two atoms, possibly periodic copies) at [`bridge_position`](@ref) with
the third atom's position `third` fixing the plane and the side. Refuses a
position closer than `min_distance` Å to any atom.
"""
function add_bridge(c::Crystal, a, b, third; element=:O, offset::Real=0.5, flip::Bool=false,
                    min_distance::Real=0.5)
    return add_atom(c, element, bridge_position(a, b, third; offset, flip); cartesian=true, min_distance)
end
