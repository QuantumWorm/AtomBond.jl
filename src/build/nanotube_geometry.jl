# Geometry of a single-wall (n, m) nanotube: a honeycomb sheet rolled into a
# cylinder. Plain Julia (LinearAlgebra only), so scripts/build_library.jl can
# include it without loading AtomBond.

"""
    nanotube_geometry(n, m, bond) -> (positions, sublattice, radius, period)

Atoms of one translational period of the (n, m) tube made from a honeycomb
sheet with nearest-neighbour distance `bond` (Å), sublattice 1 first:
Cartesian positions with
the tube axis along z through the origin, the sublattice (1 or 2) of each
atom, the tube radius and the period T along the axis (Å). The chiral vector
is C_h = n a₁ + m a₂ (a₁, a₂ at 60°), T = t₁ a₁ + t₂ a₂ with t₁ = (2m+n)/d,
t₂ = -(2n+m)/d, d = gcd(2m+n, 2n+m); one period holds 4(n²+nm+m²)/d atoms.
"""
function nanotube_geometry(n::Integer, m::Integer, bond::Real)
    (n >= 0 && m >= 0 && n + m > 0) || throw(ArgumentError("nanotube indices need n, m ≥ 0, not both 0"))
    bond > 0 || throw(ArgumentError("the bond length must be positive"))
    a = bond * sqrt(3)                                          # sheet lattice constant
    a1, a2 = [a, 0.0], [a / 2, a * sqrt(3) / 2]
    basis = [(a1 + a2) / 3, 2(a1 + a2) / 3]                     # the two sublattices
    ch = n * a1 + m * a2
    d = gcd(2m + n, 2n + m)
    t = ((2m + n) ÷ d) * a1 - ((2n + m) ÷ d) * a2
    lc, lt = norm(ch), norm(t)
    R = lc / 2π
    pos, sub = Vector{Float64}[], Int[]
    span = 2(n + m + (2m + n) ÷ d + (2n + m) ÷ d)
    for i in -span:span, j in -span:span, (k, b) in enumerate(basis)
        p = i * a1 + j * a2 + b
        u, v = dot(p, ch) / lc^2, dot(p, t) / lt^2              # fractions of C_h and T
        (-1e-9 <= u < 1 - 1e-9 && -1e-9 <= v < 1 - 1e-9) || continue
        push!(pos, [R * cos(2π * u), R * sin(2π * u), v * lt])
        push!(sub, k)
    end
    expected = 4 * (n^2 + n * m + m^2) ÷ d
    length(pos) == expected || error("($n,$m) nanotube: found $(length(pos)) atoms, expected $expected")
    order = sortperm(sub)                                       # sublattice 1 first, then 2
    return pos[order], sub[order], R, lt
end

"""Kind of (n, m) tube: armchair, zigzag or chiral, and metallic when (n - m) is a multiple of 3."""
function nanotube_kind(n::Integer, m::Integer)
    shape = n == m ? "armchair" : (n == 0 || m == 0) ? "zigzag" : "chiral"
    return shape, mod(n - m, 3) == 0
end
