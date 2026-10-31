# Lattice: three cell vectors stored as the *columns* of a 3×3 matrix (Å),
# so that   r_cart = M * r_frac   and   r_frac = M \ r_cart.

const Vec3 = SVector{3,Float64}
const Mat3 = SMatrix{3,3,Float64,9}

struct Lattice
    matrix::Mat3
    function Lattice(m::AbstractMatrix)
        size(m) == (3, 3) || throw(ArgumentError("lattice matrix must be 3×3"))
        M = Mat3(m)
        det(M) > 1e-8 || throw(ArgumentError("lattice vectors must be right-handed and non-degenerate (det = $(det(M)))"))
        new(M)
    end
end

"""
    Lattice(a, b, c)

Build a lattice from three cell vectors.
"""
Lattice(a::AbstractVector, b::AbstractVector, c::AbstractVector) = Lattice(hcat(a, b, c))

"""
    Lattice(a, b, c, α, β, γ)

Build a lattice from cell lengths (Å) and angles (degrees), using the standard
orientation: **a** along x, **b** in the xy plane.
"""
function Lattice(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real)
    cα, cβ, cγ = cosd(α), cosd(β), cosd(γ)
    sγ = sind(γ)
    cx = c * cβ
    cy = c * (cα - cβ * cγ) / sγ
    cz = sqrt(c^2 - cx^2 - cy^2)
    Lattice(Vec3(a, 0, 0), Vec3(b * cγ, b * sγ, 0), Vec3(cx, cy, cz))
end

"""Cubic lattice with edge `a`."""
cubic(a::Real) = Lattice(a, a, a, 90, 90, 90)

"""The lattice vectors a, b, c (Å)."""
cell_vectors(L::Lattice) = (L.matrix[:, 1], L.matrix[:, 2], L.matrix[:, 3])
"""Cell volume in Å³ (negative for a left-handed cell)."""
cell_volume(L::Lattice) = det(L.matrix)

"""
    lattice_parameters(L) -> (a, b, c, α, β, γ)

Lengths in Å and angles in degrees.
"""
function lattice_parameters(L::Lattice)
    a, b, c = cell_vectors(L)
    la, lb, lc = norm(a), norm(b), norm(c)
    angle(u, v) = acosd(clamp(dot(u, v) / (norm(u) * norm(v)), -1, 1))
    return (la, lb, lc, angle(b, c), angle(a, c), angle(a, b))
end

"""Cartesian position (Å) of fractional coordinates `f`."""
frac_to_cart(L::Lattice, f::AbstractVector) = L.matrix * Vec3(f)
"""Fractional coordinates of a Cartesian position `r` (Å)."""
cart_to_frac(L::Lattice, r::AbstractVector) = L.matrix \ Vec3(r)

function Base.show(io::IO, L::Lattice)
    a, b, c, α, β, γ = lattice_parameters(L)
    @printf(io, "Lattice(a=%.4f, b=%.4f, c=%.4f Å, α=%.2f, β=%.2f, γ=%.2f°)", a, b, c, α, β, γ)
end
