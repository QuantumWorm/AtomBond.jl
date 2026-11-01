# Molecules of the structure library: experimental gas-phase geometries
# (bond lengths in Å, angles in degrees) from the cited spectroscopy and
# electron-diffraction papers and books. `build_library.jl` computes the
# Cartesian coordinates from these values. Included by build_library.jl.

using LinearAlgebra

const V3 = Vector{Float64}

# --- geometry builders ----------------------------------------------------------------

"""
Position of a new atom at distance `r` from `a`, with angle new–a–b = `θ` and
dihedral new–a–b–c = `φ` (the usual Z-matrix step).
"""
function zplace(a, b, c, r, θ, φ)
    u = normalize(a - b)
    n = normalize(cross(b - c, u))
    m = cross(n, u)
    return a + r * (-cosd(θ) * u + sind(θ) * cosd(φ) * m + sind(θ) * sind(φ) * n)
end

"""
Atoms from a Z-matrix: rows (symbol,), (symbol, i, r), (symbol, i, r, j, θ) or
(symbol, i, r, j, θ, k, φ), with i, j, k earlier atoms.
"""
function zmatrix(rows)
    syms, pos = String[], V3[]
    for row in rows
        push!(syms, row[1])
        n = length(row)
        p = n == 1 ? [0.0, 0, 0] :
            n == 3 ? pos[row[2]] + [row[3], 0, 0] :
            n == 5 ? zplace(pos[row[2]], pos[row[4]], pos[row[4]] + [0, 1.0, 0], row[3], row[5], 0) :
                     zplace(pos[row[2]], pos[row[4]], pos[row[6]], row[3], row[5], row[7])
        push!(pos, p)
    end
    return syms, pos
end

"""A linear molecule: atoms along x with the given bond lengths."""
linear(syms, bonds) = (collect(syms), [[x, 0.0, 0.0] for x in cumsum([0.0; collect(bonds)])])

"""Pyramidal XY3 (C3v) with bond length r and Y–X–Y angle θ."""
function pyramid(x, y, r, θ)
    α = acos(sqrt((2cosd(θ) + 1) / 3))            # angle between a bond and the threefold axis
    ys = [r * [sin(α) * cos(2π * k / 3), sin(α) * sin(2π * k / 3), -cos(α)] for k in 0:2]
    return [x; fill(y, 3)], [[0.0, 0, 0], ys...]
end

"""Tetrahedral XY4 (Td) with bond length r."""
function tetrahedron(x, y, r)
    d = r / sqrt(3)
    return [x; fill(y, 4)], [[0.0, 0, 0], [d, d, d], [d, -d, -d], [-d, d, -d], [-d, -d, d]]
end

"""Staggered X2Y6 (D3d, ethane type) with X–X, X–Y and the X–X–Y angle β."""
function staggered(x, y, rxx, rxy, β)
    s, c = rxy * sind(β), rxy * cosd(β)
    top = [[s * cos(2π * k / 3), s * sin(2π * k / 3), rxx / 2 - c] for k in 0:2]
    bottom = [[s * cos(2π * k / 3 + π / 3), s * sin(2π * k / 3 + π / 3), -rxx / 2 + c] for k in 0:2]
    return [x, x, fill(y, 6)...], [[0.0, 0, rxx / 2], [0.0, 0, -rxx / 2], top..., bottom...]
end

"""Planar X2Y4 (D2h, ethylene type) with X=X, X–Y and the Y–X–Y angle θ."""
function planar_x2y4(x, y, rxx, rxy, θ)
    h = [[sgn * (rxx / 2 + rxy * cosd(θ / 2)), t * rxy * sind(θ / 2), 0.0] for sgn in (1, -1) for t in (1, -1)]
    return [x, x, fill(y, 4)...], [[rxx / 2, 0, 0], [-rxx / 2, 0, 0], h...]
end

"""Benzene (D6h) with C–C and C–H."""
function benzene(rcc, rch)
    ring(r) = [[r * cos(π * k / 3), r * sin(π * k / 3), 0.0] for k in 0:5]
    return [fill("C", 6); fill("H", 6)], [ring(rcc); ring(rcc + rch)]
end

"""Diborane (D2h): B–B, terminal and bridging B–H, terminal H–B–H angle."""
function diborane(rbb, rbht, rbhb, θ)
    h = sqrt(rbhb^2 - (rbb / 2)^2)
    b, ht = planar_x2y4("B", "H", rbb, rbht, θ)
    return [b; "H"; "H"], [ht; [[0.0, 0, h], [0.0, 0, -h]]]
end

"""Tetrahedral-type XY2Z2 (C2v): Y–X–Y in the xz plane, Z–X–Z in the yz plane."""
function xy2z2(x, y, z, rxy, θy, rxz, θz)
    ys = [[s * rxy * sind(θy / 2), 0.0, rxy * cosd(θy / 2)] for s in (1, -1)]
    zs = [[0.0, s * rxz * sind(θz / 2), -rxz * cosd(θz / 2)] for s in (1, -1)]
    return [x, y, y, z, z], [[0.0, 0, 0], ys..., zs...]
end

"""Trigonal planar XY3 (D3h) with bond length r, in the xy plane."""
planar_xy3(x, y, r) = ([x; fill(y, 3)], [[0.0, 0, 0], [[r * cos(2π * k / 3), r * sin(2π * k / 3), 0.0] for k in 0:2]...])

"""C3v XY3Z (methyl fluoride type): pyramidal XY3 (bond r_xy, Y–X–Y angle θ) and Z on the axis at r_xz."""
function xy3z(x, y, z, rxy, θ, rxz)
    syms, pos = pyramid(x, y, rxy, θ)
    return [syms; z], [pos; [[0.0, 0, rxz]]]
end

"""Trigonal bipyramidal XY5 (D3h) with axial and equatorial bond lengths."""
function bipyramid(x, y, rax, req)
    eq = [[req * cos(2π * k / 3), req * sin(2π * k / 3), 0.0] for k in 0:2]
    return [x; fill(y, 5)], [[0.0, 0, 0], eq..., [0.0, 0, rax], [0.0, 0, -rax]]
end

"""Octahedral XY6 (Oh) with bond length r."""
octahedron(x, y, r) = ([x; fill(y, 6)], [[0.0, 0, 0], ([s * r * (1:3 .== k) for k in 1:3 for s in (1, -1)])...])

"""
Hydrogens of a C3v methyl group on carbon `c` bonded to `m`: C–H `r`, angle
m–C–H `β`, and the first H at torsion `τ` (degrees) from the direction `ref`.
"""
function methyl_h(m, c, ref, r, β, τ)
    u = normalize(c - m)
    e1 = normalize(ref - dot(ref, u) * u)
    e2 = cross(u, e1)
    return [c + r * (-cosd(β) * u + sind(β) * (cosd(τ + 120k) * e1 + sind(τ + 120k) * e2)) for k in 0:2]
end

"""
Planar M(CH3)3 (C3h): M–C `rmc`, C–H `rch`, angle M–C–H `β`; one C–H of each
methyl lies in the MC3 plane, pointing away from the next methyl.
"""
function trimethyl(x, rmc, rch, β)
    cs = [rmc * [cos(2π * k / 3), sin(2π * k / 3), 0.0] for k in 0:2]
    hs = [h for k in 1:3 for h in methyl_h([0.0, 0, 0], cs[k], cs[mod1(k + 1, 3)], rch, β, 180.0)]
    return [x; fill("C", 3); fill("H", 9)], [[[0.0, 0, 0]]; cs; hs]
end

"""Linear M(CH3)2 with staggered methyl groups (D3d)."""
function dimethyl(x, rmc, rch, β)
    c1, c2 = [0.0, 0, rmc], [0.0, 0, -rmc]
    hs = [methyl_h([0.0, 0, 0], c1, [1.0, 0, 0], rch, β, 0.0); methyl_h([0.0, 0, 0], c2, [1.0, 0, 0], rch, β, 60.0)]
    return [x, "C", "C", fill("H", 6)...], [[[0.0, 0, 0], c1, c2]; hs]
end

"""Metal carbonyl M(CO)n, n = 4 (tetrahedral) or 6 (octahedral), with M–C and C–O."""
function carbonyl(x, n, rmc, rco)
    dirs = n == 4 ? normalize.([[1.0, 1, 1], [1.0, -1, -1], [-1.0, 1, -1], [-1.0, -1, 1]]) :
                    [Float64.(s * (1:3 .== k)) for k in 1:3 for s in (1, -1)]
    return [x; fill("C", n); fill("O", n)], [[[0.0, 0, 0]]; [rmc * d for d in dirs]; [(rmc + rco) * d for d in dirs]]
end

"""
Methyl-bridged M2(CH3)6 (trimethylaluminium dimer): M–M, terminal and bridging
M–C, terminal C–M–C angle, C–H and M–C–H (also used, as the angle to the
outward axis, for the bridging methyls). Terminal methyls staggered with
respect to the other bonds of their metal; bridging methyls with one C–H in
the plane through both carbons perpendicular to M–M.
"""
function bridged_trimethyl_dimer(x, rmm, rt, rb, θt, rch, β)
    ms = [[rmm / 2, 0, 0], [-rmm / 2, 0, 0]]
    yb = sqrt(rb^2 - (rmm / 2)^2)
    cb = [[0.0, yb, 0], [0.0, -yb, 0]]
    ct = [[sg * (rmm / 2 + rt * cosd(θt / 2)), 0.0, t * rt * sind(θt / 2)] for sg in (1, -1) for t in (1, -1)]
    other = [ct[2], ct[1], ct[4], ct[3]]                      # the other terminal carbon on the same metal
    ht = [h for k in 1:4 for h in methyl_h(ms[(k + 1) ÷ 2], ct[k], other[k] - ms[(k + 1) ÷ 2], rch, β, 180.0)]
    hb = [h for c in cb for h in methyl_h(-c, c, [0.0, 0, 1], rch, β, 0.0)]  # axis outward from the M2 centre
    return [x, x, fill("C", 6)..., fill("H", 18)...], [ms; cb; ct; ht; hb]
end

"""
Halogen-bridged M2X6 (Al2Cl6 type, D2h): bridging and terminal M–X, X(b)–M–X(b)
and X(t)–M–X(t) angles. The M2X2 ring lies in the xy plane.
"""
function bridged_halide_dimer(x, y, rb, φb, rt, θt)
    a = rb * cosd(φb / 2)
    ms = [[a, 0.0, 0], [-a, 0.0, 0]]
    xb = [[0.0, rb * sind(φb / 2), 0], [0.0, -rb * sind(φb / 2), 0]]
    xt = [[sg * (a + rt * cosd(θt / 2)), 0.0, t * rt * sind(θt / 2)] for sg in (1, -1) for t in (1, -1)]
    return [x, x, fill(y, 6)...], [ms; xb; xt]
end

"""Atoms at given Cartesian coordinates (Å), e.g. a paper's substitution coordinates."""
cartesian(rows) = ([r[1] for r in rows], [Float64[r[2], r[3], r[4]] for r in rows])

"""
Three-membered ring X–C2 with CH2 groups (oxirane type, C2v): C–C, C–X, C–H,
H–C–H, and θ, the angle between each H2C plane and the C–C bond (the H atoms
bent away from X).
"""
function ring_ch2(x, rcc, rcx, rch, θhch, θ)
    yx = sqrt(rcx^2 - (rcc / 2)^2)
    cs = [[rcc / 2, 0.0, 0], [-rcc / 2, 0.0, 0]]
    hs = V3[]
    for (sg, c) in zip((1, -1), cs)
        b = [sg * cosd(θ), -sind(θ), 0.0]                      # bisector, away from the ring and from X
        for t in (1, -1)
            push!(hs, c + rch * (cosd(θhch / 2) * b + t * sind(θhch / 2) * [0.0, 0, 1]))
        end
    end
    return ["C", "C", x, "H", "H", "H", "H"], [cs; [[0.0, yx, 0]]; hs]
end

"""C3v methyl hydrogens' angle to their axis from H–C–H (the same relation as for a pyramid)."""
methyl_axis_angle(θhch) = acosd(sqrt((2cosd(θhch) + 1) / 3))

"""
C3v XY3Z from the Z–X–Y angle (chloromethane, chloroform): X–Y, X–Z and the
angle between the X–Z bond and each X–Y bond.
"""
function xy3z_axial(x, y, z, rxy, θzxy, rxz)
    α = 180 - θzxy                                             # Y bonds from the axis pointing away from Z
    ys = [rxy * [sind(α) * cos(2π * k / 3), sind(α) * sin(2π * k / 3), -cosd(α)] for k in 0:2]
    return [x; fill(y, 3); z], [[[0.0, 0, 0]]; ys; [[0.0, 0, rxz]]]
end

"""
In-plane (z = 0) atom at distance `r` from `a` with angle new–a–b = `θ`; of the
two such points, the one farther from `away`.
"""
function place2d(a, b, r, θ, away)
    u = normalize(b - a)
    n = [-u[2], u[1], 0.0]
    p1, p2 = a + r * (cosd(θ) * u + sind(θ) * n), a + r * (cosd(θ) * u - sind(θ) * n)
    return norm(p1 - away) > norm(p2 - away) ? p1 : p2
end

mirror_y(p) = [p[1], -p[2], p[3]]

"""
Planar five-membered ring X1 C2 C3 C4 C5 with H on each carbon (furan,
thiophene; C2v): X–C2, C2–C3, C2–X–C5 angle, X–C2–C3 angle, C–H lengths and
the angles X–C2–H and C2–C3–H.
"""
function ring5(x, rxc, rcc, θcxc, θxcc, rch2, θxch, rch3, θcch)
    X = [0.0, 0, 0]
    c2 = rxc * [cosd(θcxc / 2), sind(θcxc / 2), 0]
    c3 = place2d(c2, X, rcc, θxcc, [-10.0, 0, 0])
    mid = (c2 + c3 + mirror_y(c2) + mirror_y(c3) + X) / 5
    h2 = place2d(c2, X, rch2, θxch, mid)
    h3 = place2d(c3, c2, rch3, θcch, mid)
    return [x, "C", "C", "C", "C", "H", "H", "H", "H"],
           [X, c2, c3, mirror_y(c3), mirror_y(c2), h2, h3, mirror_y(h3), mirror_y(h2)]
end

"""
Pyridine (C2v) from its ring bonds and angles and the three C–H bonds with
their angles C3–C2–H, C4–C3–H and C3–C4–H.
"""
function pyridine(rnc, r23, r34, θcnc, θncc, θccc, rch, θh)
    N = [0.0, 0, 0]
    c2 = rnc * [cosd(θcnc / 2), sind(θcnc / 2), 0]
    c3 = place2d(c2, N, r23, θncc, [-10.0, 0, 0])
    c4 = place2d(c3, c2, r34, θccc, [c3[1], 100.0, 0])                # toward the mirror line (y = 0)
    mid = (N + c4) / 2
    h2 = place2d(c2, c3, rch[1], θh[1], mid)
    h3 = place2d(c3, c4, rch[2], θh[2], mid)
    h4 = place2d(c4, c3, rch[3], θh[3], mid)
    return ["N", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H"],
           [N, c2, c3, c4, mirror_y(c3), mirror_y(c2), h2, h3, h4, mirror_y(h3), mirror_y(h2)]
end

"""Acetonitrile-type H3C–C≡N: C–C, C≡N, C–H and the C–C–H angle."""
function methyl_nitrile(rcc, rcn, rch, θcch)
    c1, c2, n = [0.0, 0, 0], [rcc, 0.0, 0], [rcc + rcn, 0.0, 0]
    return ["C", "C", "N", "H", "H", "H"], [c1, c2, n, methyl_h(c2, c1, [0.0, 1, 0], rch, θcch, 0.0)...]
end

"""Propyne-type H3C–C≡C–H: C–C, C≡C, methyl C–H, H–C–C angle and acetylenic C–H."""
function propyne(rcc, rcc3, rch3, θhcc, rch)
    c1, c2 = [0.0, 0, 0], [rcc, 0.0, 0]
    c3, h = c2 + [rcc3, 0, 0], c2 + [rcc3 + rch, 0, 0]
    return ["C", "C", "C", "H", "H", "H", "H"], [c1, c2, c3, h, methyl_h(c2, c1, [0.0, 1, 0], rch3, θhcc, 0.0)...]
end

"""Allene H2C=C=CH2 (D2d): C=C, C–H, H–C=C; the two CH2 planes perpendicular."""
function allene(rcc, rch, θ)
    hx = rcc + rch * cosd(180 - θ)
    hs = [[hx, s * rch * sind(θ), 0.0] for s in (1, -1)]
    return ["C", "C", "C", "H", "H", "H", "H"],
           [[0.0, 0, 0], [rcc, 0.0, 0], [-rcc, 0.0, 0], hs..., [[-hx, 0.0, s * rch * sind(θ)] for s in (1, -1)]...]
end

"""Cyclopropane (D3h): C–C, C–H, H–C–H; each HCH plane perpendicular to the ring."""
function cyclopropane(rcc, rch, θhch)
    R = rcc / sqrt(3)
    cs = [R * [cos(2π * k / 3), sin(2π * k / 3), 0.0] for k in 0:2]
    hs = [c + rch * (cosd(θhch / 2) * normalize(c) + t * sind(θhch / 2) * [0.0, 0, 1]) for c in cs for t in (1, -1)]
    return [fill("C", 3); fill("H", 6)], [cs; hs]
end

"""
Puckered cyclobutane (D2d): C–C, C–H (axial = equatorial), H–C–H, the ring
dihedral θ (angle between the planes C2C1C4 and C2C3C4) and the CH2 rocking
angle β (between the HCH and CCC bisectors), the axial C–H bonds in 1,3
positions tilted towards each other.
"""
function cyclobutane(rcc, rch, θhch, θ, β)
    a = rcc / sqrt(2 + tand(θ / 2)^2)
    h = a * tand(θ / 2) / 2
    c1 = [a, 0.0, h]
    u = normalize([a, 0.0, 2h])                                # CCC bisector at C1, outward
    rot(v, ang) = [cosd(ang) * v[1] + sind(ang) * v[3], v[2], -sind(ang) * v[1] + cosd(ang) * v[3]]
    function hydrogens(sgn)
        b = rot(u, sgn * β)
        n = normalize(cross([0.0, 1, 0], b))                   # in the xz plane, perpendicular to b
        n[3] < 0 && (n = -n)
        return c1 + rch * (cosd(θhch / 2) * b + sind(θhch / 2) * n), c1 + rch * (cosd(θhch / 2) * b - sind(θhch / 2) * n)
    end
    hax, heq = argmin(p -> p[1][1], [hydrogens(1), hydrogens(-1)])   # axial H leans towards C3 (−x)
    s4(p) = [p[2], -p[1], -p[3]]                               # the S4 operation maps C1 → C4 → C3 → C2
    cs, hs = [c1], [hax, heq]
    for _ in 1:3
        push!(cs, s4(cs[end])); append!(hs, s4.(hs[end-1:end]))
    end
    return [fill("C", 4); fill("H", 8)], [cs; hs]
end

"""
Buckminsterfullerene C60 (Ih, truncated icosahedron) from its two bonds: `r56`
(in a pentagon, fusing a 5- and a 6-ring) and `r66` (fusing two 6-rings). Each
atom lies on an edge of an icosahedron with edge r66 + 2 r56, at r56 from a vertex.
"""
function fullerene_c60(r56, r66)
    φ = (1 + sqrt(5)) / 2
    ico = [Float64[0, s1, s2 * φ] for s1 in (1, -1) for s2 in (1, -1)]
    ico = [ico; [circshift(v, 1) for v in ico]; [circshift(v, 2) for v in ico]]   # 12 vertices, edge 2
    L = r66 + 2r56
    ico = ico .* (L / 2)
    pos = V3[]
    for i in eachindex(ico), j in eachindex(ico)
        i != j && isapprox(norm(ico[i] - ico[j]), L; rtol=1e-9) && push!(pos, ico[i] + r56 / L * (ico[j] - ico[i]))
    end
    return fill("C", length(pos)), pos
end

"""
Sulfuric acid (C2): S–O(H), S=O, O–H, the angles O1–S–O1′, O2–S–O2′ and
H–O1–S, the angle between the O1SO1′ and O2SO2′ planes (both bisected by the
C2 axis) and the torsion O1′–S–O1–H.
"""
function sulfuric_acid(rso1, rso2, roh, θ11, θ22, planes, θhos, τ)
    pair(r, θ, φ, up) = (s = r * sind(θ / 2); z = up * r * cosd(θ / 2);
                         ([s * cosd(φ), s * sind(φ), z], [-s * cosd(φ), -s * sind(φ), z]))
    o1, o1p = pair(rso1, θ11, 0.0, 1)
    o2, o2p = pair(rso2, θ22, planes, -1)
    S = [0.0, 0, 0]
    h = zplace(o1, S, o1p, roh, θhos, τ)
    hp = [-h[1], -h[2], h[3]]                                   # C2 about z
    return ["S", "O", "O", "O", "O", "H", "H"], [S, o1, o1p, o2, o2p, h, hp]
end

"""Cs XYZ2 (thionyl chloride type): X–Y, X–Z, the Y–X–Z angle (both Z) and the Z–X–Z angle."""
function xyz2_pyramid(x, y, z, rxy, rxz, θyxz, θzxz)
    φ = acosd((cosd(θzxz) - cosd(θyxz)^2) / sind(θyxz)^2)        # azimuth between the two Z about X–Y
    zs = [rxz * [sind(θyxz) * cosd(s * φ / 2), sind(θyxz) * sind(s * φ / 2), cosd(θyxz)] for s in (1, -1)]
    return [x, y, z, z], [[0.0, 0, 0], [0.0, 0, rxy], zs...]
end

"""Tetrahedral X4 (P4, Td) with edge (bond) length r."""
function tetrahedral_x4(x, r)
    d = r / (2sqrt(2))
    return fill(x, 4), [[d, d, d], [d, -d, -d], [-d, d, -d], [-d, -d, d]]
end

"""
Crown ring X_n (S8, D_{n/2 d}): atoms alternately up and down on two circles,
from the bond length r and the bond angle θ (the torsion follows).
"""
function crown_ring(x, n, r, θ)
    β = 2π / n                                                  # azimuth step
    # Bond length and angle fix the radius R and the half-height h:
    #   r² = 2R²(1 - cos β) + 4h²,  angle k-1, k, k+1 (neighbours at the same height)
    #   cos θ = (2R²(cos β - cos 2β)... solved numerically below.
    pos(R, h) = [[R * cos(β * k), R * sin(β * k), iseven(k) ? h : -h] for k in 0:n-1]
    ang(R, h) = (p = pos(R, h); u = p[n] - p[1]; v = p[2] - p[1]; acosd(dot(u, v) / (norm(u) * norm(v))))
    # For a given R, h follows from r; find R by bisection on the angle.
    hof(R) = sqrt(max(r^2 - 2R^2 * (1 - cos(β)), 0.0)) / 2
    lo, hi = 0.5r, r / sqrt(2(1 - cos(β)))                     # h > 0 needs R below hi
    for _ in 1:200
        mid = (lo + hi) / 2
        ang(mid, hof(mid)) < θ ? (lo = mid) : (hi = mid)
    end
    R = (lo + hi) / 2
    return fill(x, n), pos(R, hof(R))
end

"""
Methylsilane-type H3C–XH3 (C3v, staggered): C–X, C–H, X–H and the H–C–H and
H–X–H angles (both groups C3v about the C–X axis).
"""
function methyl_xh3(x, rcx, rch, rxh, θhch, θhxh)
    c, xx = [0.0, 0, 0], [0.0, 0, rcx]
    hc = methyl_h(xx, c, [1.0, 0, 0], rch, 180 - methyl_axis_angle(θhch), 0.0)
    hx = methyl_h(c, xx, [1.0, 0, 0], rxh, 180 - methyl_axis_angle(θhxh), 60.0)
    return ["C", x, "H", "H", "H", "H", "H", "H"], [[c, xx]; hc; hx]
end

"""
Disiloxane-type H3X–O–XH3 (C2v): X–O, X–H, X–O–X, H–X–O (XH3 groups C3v
about X–O), with one H of each group in the XOX plane, towards the other X.
"""
function disiloxane(x, rxo, rxh, θxox, θhxo)
    o = [0.0, 0, 0]
    x1 = rxo * [sind(θxox / 2), 0, cosd(θxox / 2)]
    x2 = rxo * [-sind(θxox / 2), 0, cosd(θxox / 2)]
    h1 = methyl_h(o, x1, x2 - x1, rxh, θhxo, 0.0)
    h2 = methyl_h(o, x2, x1 - x2, rxh, θhxo, 0.0)
    return ["O", x, x, fill("H", 6)...], [[o, x1, x2]; h1; h2]
end

"""
Trisilane-type H3X–XH2–XH3 (C2v, terminal groups staggered): X–X, X–H
(terminal and central), X–X–X, X(2)–X(1)–H and H–X(2)–H; one terminal H in
the XXX plane, anti to the far X.
"""
function trisilane(x, rxx, rxh1, rxh2, θxxx, θxxh, θhxh)
    x2 = [0.0, 0, 0]
    x1 = rxx * [sind(θxxx / 2), 0, -cosd(θxxx / 2)]
    x3 = rxx * [-sind(θxxx / 2), 0, -cosd(θxxx / 2)]
    hc = [rxh2 * [0, s * sind(θhxh / 2), cosd(θhxh / 2)] for s in (1, -1)]
    h1 = methyl_h(x2, x1, x1 - x3, rxh1, θxxh, 0.0)                # first H anti to X(3)
    h3 = methyl_h(x2, x3, x3 - x1, rxh1, θxxh, 0.0)
    return [x, x, x, fill("H", 8)...], [[x1, x2, x3]; hc; h1; h3]
end

# --- references --------------------------------------------------------------------------

const HUBER_HERZBERG = ("K. P. Huber, G. Herzberg, Molecular Spectra and Molecular Structure IV. Constants of " *
                        "Diatomic Molecules, Van Nostrand Reinhold (1979)", "10.1007/978-1-4757-0961-2")
const HERZBERG_1966 = ("G. Herzberg, Molecular Spectra and Molecular Structure III. Electronic Spectra and " *
                       "Electronic Structure of Polyatomic Molecules, Van Nostrand (1966)", "")
const KUCHITSU = ("K. Kuchitsu (ed.), Structure of Free Polyatomic Molecules - Basic Data, Springer (1998)",
                  "10.1007/978-3-642-45748-7")
const GURVICH = ("L. V. Gurvich, I. V. Veyts, C. B. Alcock, Thermodynamic Properties of Individual Substances, " *
                 "4th ed., Hemisphere (1989)", "")

# Formamide: N–H' (trans to the formyl H) at dihedral H'–N–C–O = 7°, and N–H''
# at H''–N–C–H = 12°, both to the same side of the O=C–N plane (Costain & Dowling).
const FORMAMIDE_TWIST = (7.0, -12.0)

"""Methyl formate (cis), from Curl 1959: planar heavy atoms, methyl staggering the C–O bond."""
function methyl_formate()
    syms, pos = zmatrix([("C",), ("O", 1, 1.200), ("O", 1, 1.334, 2, 125.867), ("H", 1, 1.101, 3, 109.3, 2, 180.0),
                         ("C", 3, 1.437, 1, 114.783, 2, 0.0)])
    hs = methyl_h(pos[3], pos[5], pos[1] - pos[3], 1.086, 180 - methyl_axis_angle(110.667), 180.0)
    return [syms; fill("H", 3)], [pos; hs]
end

"""Pentagonal bipyramidal XY7 (D5h) with axial and equatorial bond lengths."""
function pentagonal_bipyramid(x, y, rax, req)
    eq = [[req * cos(2π * k / 5), req * sin(2π * k / 5), 0.0] for k in 0:4]
    return [x; fill(y, 7)], [[0.0, 0, 0], eq..., [0.0, 0, rax], [0.0, 0, -rax]]
end

"""
Planar molecule built on the bond A–B (A at the origin, B along x):
substituents (symbol, 1 = on A or 2 = on B, bond length, angle to the A–B
bond, side +1 or -1 of the bond) in the xy plane.
"""
function planar_ab(a, b, rab, subs)
    pa, pb = [0.0, 0, 0], [rab, 0.0, 0]
    syms, pos = [a, b], V3[pa, pb]
    for (s, at, r, θ, side) in subs
        o, other = at == 1 ? (pa, pb) : (pb, pa)
        push!(syms, s)
        push!(pos, place2d(o, other, r, θ, [0.0, -side * 100.0, 0]))
    end
    return syms, pos
end

"""
Cyclopropene (C2v): C1–C2, C2=C2′, methylene and ethylenic C–H, the angle
C2′=C2–H and H–C1–H; the CH2 plane perpendicular to the ring.
"""
function cyclopropene(r12, r22, rch2, rch1, θh, θhch)
    c2s = [[r22 / 2, 0.0, 0], [-r22 / 2, 0.0, 0]]
    c1 = [0.0, sqrt(r12^2 - (r22 / 2)^2), 0]
    hs = [place2d(c2s[1], c2s[2], rch1, θh, c1), place2d(c2s[2], c2s[1], rch1, θh, c1)]
    hm = [c1 + rch2 * (cosd(θhch / 2) * [0.0, 1, 0] + t * sind(θhch / 2) * [0.0, 0, 1]) for t in (1, -1)]
    return ["C", "C", "C", "H", "H", "H", "H"], [[c1]; c2s; hs; hm]
end

"""
Cs XZY2W (CHClF2 type): Z along +z at r_xz, two Y at angle Z–X–Y with
Y–X–Y between them, and W in the mirror plane at angle W–X–Z on the side
away from the Y atoms.
"""
function xzy2w(x, z, y, w, rxz, rxy, rxw, θzxy, θyxy, θwxz)
    φ = acos((cosd(θyxy) - cosd(θzxy)^2) / sind(θzxy)^2)       # azimuth between the two Y bonds
    ys = [rxy * [sind(θzxy) * cos(s * φ / 2), sind(θzxy) * sin(s * φ / 2), cosd(θzxy)] for s in (1, -1)]
    return [x, z, y, y, w], [[0.0, 0, 0], [0.0, 0, rxz], ys..., rxw * [-sind(θwxz), 0.0, cosd(θwxz)]]
end

"""Angle in degrees from degrees and arc minutes."""
dm(d, m) = d + m / 60

"""Note for singlet carbenes and silylenes (closed shell)."""
const SINGLET_NOTE = "Ground state 1A1 (singlet, closed shell)."

const TYLER_SHERIDAN = ("J. K. Tyler, J. Sheridan, Structural studies of linear molecules by microwave spectroscopy, " *
                        "Trans. Faraday Soc. 59, 2661 (1963); values read in Table 5 of the authors' report of the " *
                        "same title (USAF contract AF61(052)-241, Feb. 1963, DTIC AD0409439, " *
                        "doi:10.21236/ad0409439)", "10.1039/tf9635902661")
const TYLER_SHERIDAN_NOTE = "Substitution structure r_s (authors' estimate: within 0.002 Å of r_e)."

spin_note(state) = "The ground state is a spin $state (open shell): use a spin-polarized calculation."
odd_electron_note(extra="") = "A radical with an odd number of electrons: the ground state is a spin doublet (open " *
                              "shell), use a spin-polarized calculation." * (isempty(extra) ? "" : " " * extra)

# name, formula, category, point group, geometry (symbols, positions), references
# (text, DOI), optional note. The values are those of the references.
const MOLECULES = [
    # Diatomic molecules (equilibrium bond lengths r_e)
    (name="Hydrogen", formula="H2", category="Inorganic: elements", pg="D∞h",
     geometry=linear(("H", "H"), (0.7414,)), refs=[HUBER_HERZBERG]),
    (name="Nitrogen", formula="N2", category="Inorganic: elements", pg="D∞h",
     geometry=linear(("N", "N"), (1.0977,)), refs=[HUBER_HERZBERG]),
    (name="Oxygen", formula="O2", category="Inorganic: elements", pg="D∞h",
     geometry=linear(("O", "O"), (1.2075,)), refs=[HUBER_HERZBERG], note=spin_note("triplet")),
    (name="Carbon monoxide", formula="CO", category="Inorganic: oxides and sulfides", pg="C∞v",
     geometry=linear(("C", "O"), (1.1283,)), refs=[HUBER_HERZBERG]),
    (name="Nitric oxide", formula="NO", category="Inorganic: oxides and sulfides", pg="C∞v",
     geometry=linear(("N", "O"), (1.1508,)), refs=[HUBER_HERZBERG], note=spin_note("doublet")),
    (name="Hydrogen fluoride", formula="HF", category="Inorganic: hydrides", use="etch gas", pg="C∞v",
     geometry=linear(("H", "F"), (0.9168,)),
     refs=[("R. J. Le Roy, Improved parameterization for combined isotopomer analysis of diatomic spectra and its " *
            "application to HF and DF, J. Mol. Spectrosc. 194, 189 (1999)", "10.1006/jmsp.1998.7786")]),
    (name="Hydrogen chloride", formula="HCl", category="Inorganic: hydrides", use="etch gas", pg="C∞v",
     geometry=linear(("H", "Cl"), (1.2746,)), refs=[HUBER_HERZBERG]),
    (name="Fluorine", formula="F2", category="Inorganic: elements", use="etch gas", pg="D∞h",
     geometry=linear(("F", "F"), (1.41426,)),
     refs=[("T. R. Gilson, M. R. Purnell, Improved rotational constants for fluorine 19F2 using coherent anti-Stokes " *
            "Raman spectroscopy of the vibrational Q branch, J. Raman Spectrosc. 12, 222 (1982)", "10.1002/jrs.1250120304")]),
    (name="Chlorine", formula="Cl2", category="Inorganic: elements", use="etch gas", pg="D∞h",
     geometry=linear(("Cl", "Cl"), (1.9882,)),
     refs=[("H. G. M. Edwards, E. A. M. Good, D. A. Long, Pure rotational Raman spectra of the chlorine species, 35Cl2 " *
            "and 35Cl37Cl, J. Chem. Soc., Faraday Trans. 2 72, 927 (1976)", "10.1039/F29767200927")],
     note="r_e of 35Cl2."),
    (name="Bromine", formula="Br2", category="Inorganic: elements", pg="D∞h",
     geometry=linear(("Br", "Br"), (2.28102,)),
     refs=[("C. Focsa, H. Li, P. F. Bernath, Characterization of the ground state of Br2 by laser-induced fluorescence " *
            "Fourier transform spectroscopy of the B3Π0+u–X1Σ+g system, J. Mol. Spectrosc. 200, 104 (2000)",
            "10.1006/jmsp.1999.8039")]),
    (name="Hydrogen bromide", formula="HBr", category="Inorganic: hydrides", use="etch gas", pg="C∞v",
     geometry=linear(("H", "Br"), (1.4146,)),
     refs=[("G. Jones, W. Gordy, Submillimeter-wave spectra of HCl and HBr, Phys. Rev. 136, A1229 (1964)",
            "10.1103/PhysRev.136.A1229")]),
    # Small molecules
    (name="Water", formula="H2O", category="Inorganic: hydrides", use="ALD/CVD co-reactant (oxidant, nitrogen source)", pg="C2v",
     geometry=zmatrix([("O",), ("H", 1, 0.9578), ("H", 1, 0.9578, 2, 104.48)]),
     refs=[("A. R. Hoy, P. R. Bunker, A precise solution of the rotation bending Schrödinger equation for a " *
            "triatomic molecule with application to the water molecule, J. Mol. Spectrosc. 74, 1 (1979)",
            "10.1016/0022-2852(79)90019-5"), HERZBERG_1966]),
    (name="Ammonia", formula="NH3", category="Inorganic: hydrides", use="ALD/CVD co-reactant (oxidant, nitrogen source)", pg="C3v",
     geometry=pyramid("N", "H", 1.0124, 106.67), refs=[HERZBERG_1966]),
    (name="Carbon dioxide", formula="CO2", category="Inorganic: oxides and sulfides", pg="D∞h",
     geometry=linear(("O", "C", "O"), (1.1621, 1.1621)), refs=[HERZBERG_1966]),
    (name="Hydrogen sulfide", formula="H2S", category="Inorganic: hydrides", pg="C2v",
     geometry=zmatrix([("S",), ("H", 1, 1.3356), ("H", 1, 1.3356, 2, 92.11)]),
     refs=[("R. L. Cook, F. C. De Lucia, P. Helminger, Molecular force field and structure of hydrogen sulfide: " *
            "recent microwave results, J. Mol. Struct. 28, 237 (1975)", "10.1016/0022-2860(75)80094-9")]),
    (name="Nitrous oxide", formula="N2O", category="Inorganic: oxides and sulfides", use="ALD/CVD co-reactant (oxidant, nitrogen source)", pg="C∞v",
     geometry=linear(("N", "N", "O"), (1.1282, 1.1842)), refs=[HERZBERG_1966]),
    (name="Sulfur dioxide", formula="SO2", category="Inorganic: oxides and sulfides", pg="C2v",
     geometry=zmatrix([("S",), ("O", 1, 1.4321), ("O", 1, 1.4321, 2, 119.5)]), refs=[HERZBERG_1966]),
    (name="Ozone", formula="O3", category="Inorganic: elements", use="ALD/CVD co-reactant (oxidant, nitrogen source)", pg="C2v",
     geometry=zmatrix([("O",), ("O", 1, 1.278), ("O", 1, 1.278, 2, 116.8)]), refs=[HERZBERG_1966]),
    (name="Hydrogen peroxide", formula="H2O2", category="Inorganic: hydrides", pg="C2",
     geometry=zmatrix([("O",), ("O", 1, 1.475), ("H", 1, 0.950, 2, 94.8), ("H", 2, 0.950, 1, 94.8, 3, 119.8)]),
     refs=[("R. L. Redington, W. B. Olson, P. C. Cross, Studies of hydrogen peroxide: the infrared spectrum and " *
            "the internal rotation problem, J. Chem. Phys. 36, 1311 (1962)", "10.1063/1.1732733")]),
    (name="Hydrogen cyanide", formula="HCN", category="Organic: nitrogen and sulfur compounds", pg="C∞v",
     geometry=linear(("H", "C", "N"), (1.064, 1.156)), refs=[HERZBERG_1966]),
    (name="Carbonyl sulfide", formula="OCS", category="Inorganic: oxides and sulfides", pg="C∞v",
     geometry=linear(("O", "C", "S"), (1.155386, 1.562021)),
     refs=[("A. Foord, J. G. Smith, D. H. Whiffen, An anharmonic force field for carbonyl sulphide, Mol. Phys. 29, " *
            "1685 (1975)", "10.1080/00268977500101491")]),
    # Organic molecules
    (name="Methane", formula="CH4", category="Organic: hydrocarbons", pg="Td",
     geometry=tetrahedron("C", "H", 1.087),
     refs=[("E. Hirota, Anharmonic potential function and equilibrium structure of methane, J. Mol. Spectrosc. " *
            "77, 213 (1979)", "10.1016/0022-2852(79)90103-6"), HERZBERG_1966]),
    (name="Ethane", formula="C2H6", category="Organic: hydrocarbons", pg="D3d",
     geometry=staggered("C", "H", 1.536, 1.091, 110.91), refs=[HERZBERG_1966, KUCHITSU]),
    (name="Ethylene", formula="C2H4", category="Organic: hydrocarbons", pg="D2h",
     geometry=planar_x2y4("C", "H", 1.339, 1.086, 117.6), refs=[HERZBERG_1966]),
    (name="Acetylene", formula="C2H2", category="Organic: hydrocarbons", pg="D∞h",
     geometry=linear(("H", "C", "C", "H"), (1.0631, 1.2026, 1.0631)), refs=[HERZBERG_1966, KUCHITSU]),
    (name="Benzene", formula="C6H6", category="Organic: hydrocarbons", pg="D6h",
     geometry=benzene(1.397, 1.084), refs=[HERZBERG_1966]),
    (name="Methanol", formula="CH3OH", category="Organic: alcohols and ethers", use="solvent", pg="Cs",
     geometry=zmatrix([("C",), ("O", 1, 1.427), ("H", 2, 0.9561, 1, 108.87),
                       ("H", 1, 1.096, 2, 109.91, 3, 180.0), ("H", 1, 1.096, 2, 109.91, 3, 60.0),
                       ("H", 1, 1.096, 2, 109.91, 3, -60.0)]),
     refs=[("P. Venkateswarlu, W. Gordy, Methyl alcohol. II. Molecular structure, J. Chem. Phys. 23, 1200 (1955)",
            "10.1063/1.1742240"),
           ("E. Herbst, J. K. Messer, F. C. De Lucia, A new analysis and additional measurements of the millimeter " *
            "and submillimeter spectrum of methanol, J. Mol. Spectrosc. 108, 42 (1984)", "10.1016/0022-2852(84)90285-6")]),
    (name="Ethanol", formula="C2H5OH", category="Organic: alcohols and ethers", use="solvent", pg="Cs",
     geometry=zmatrix([("C",), ("C", 1, 1.512), ("O", 2, 1.4311, 1, 107.80), ("H", 3, 0.971, 2, 105.40, 1, 180.0),
                       ("H", 1, 1.098, 2, 107.81, 3, 180.0), ("H", 1, 1.088, 2, 107.82, 3, -58.90),
                       ("H", 1, 1.088, 2, 107.82, 3, 58.90), ("H", 2, 1.0859, 3, 107.83, 1, -121.10),
                       ("H", 2, 1.0859, 3, 107.83, 1, 121.10)]),
     refs=[("R. K. Kakar, C. R. Quade, Microwave rotational spectrum and internal rotation in gauche ethyl alcohol, " *
            "J. Chem. Phys. 72, 4300 (1980)", "10.1063/1.439723")],
     note="The trans (anti) conformer, a model geometry: check the parameters against the literature before " *
          "relying on them (the cited microwave study is of the gauche conformer)."),
    (name="Formaldehyde", formula="H2CO", category="Organic: carbonyl compounds", pg="C2v",
     geometry=zmatrix([("C",), ("O", 1, 1.205), ("H", 1, 1.111, 2, 121.93), ("H", 1, 1.111, 2, 121.93, 3, 180.0)]),
     refs=[HERZBERG_1966, GURVICH]),
    (name="Formic acid", formula="HCOOH", category="Organic: carbonyl compounds", pg="Cs",
     geometry=zmatrix([("C",), ("O", 1, 1.202), ("O", 1, 1.343, 2, 124.90), ("H", 1, 1.097, 3, 111.0, 2, 180.0),
                       ("H", 3, 0.972, 1, 106.3, 2, 0.0)]),
     refs=[HERZBERG_1966]),
    # Process gases and precursors (CVD, ALD, etching)
    (name="Silane", formula="SiH4", category="Inorganic: hydrides", use="deposition precursor", pg="Td",
     geometry=tetrahedron("Si", "H", 1.4799),
     refs=[("D. R. J. Boyd, Infrared spectrum of trideuterosilane and the structure of the silane molecule, " *
            "J. Chem. Phys. 23, 922 (1955)", "10.1063/1.1742148")]),
    (name="Disilane", formula="Si2H6", category="Inorganic: hydrides", use="deposition precursor", pg="D3d",
     geometry=staggered("Si", "H", 2.32, 1.47, 109.40),
     refs=[("G. W. Bethke, M. K. Wilson, Vibrational spectrum of disilane, J. Chem. Phys. 26, 1107 (1957)",
            "10.1063/1.1743480")],
     note="An approximate structure from the vibrational analysis (nearly tetrahedral angles)."),
    (name="Phosphine", formula="PH3", category="Inorganic: hydrides", use="dopant source", pg="C3v",
     geometry=pyramid("P", "H", 1.421, 93.30), refs=[HERZBERG_1966]),
    (name="Arsine", formula="AsH3", category="Inorganic: hydrides", use="dopant source", pg="C3v",
     geometry=pyramid("As", "H", 1.5108, 92.08),
     refs=[("G. S. Blevins, A. W. Jache, W. Gordy, Millimeter wave spectra of AsH3 and AsD3, Phys. Rev. 97, 684 " *
            "(1955)", "10.1103/PhysRev.97.684")]),
    (name="Diborane", formula="B2H6", category="Inorganic: hydrides", use="dopant source", pg="D2h",
     geometry=diborane(1.76, 1.20, 1.317, 122.0), refs=[GURVICH]),
    (name="Dichlorosilane", formula="SiH2Cl2", category="Inorganic: nonmetal halides", use="deposition precursor", pg="C2v",
     geometry=xy2z2("Si", "H", "Cl", 1.480, 111.30, 2.033, 109.72),
     refs=[("R. W. Davis, M. C. L. Gerry, The microwave spectrum, structure, chlorine nuclear quadrupole coupling " *
            "constants, dipole moment, and centrifugal distortion constants of dichlorosilane, J. Mol. Spectrosc. " *
            "60, 117 (1976)", "10.1016/0022-2852(76)90120-X")]),
    (name="Silicon tetrachloride", formula="SiCl4", category="Inorganic: nonmetal halides", use="deposition precursor", pg="Td",
     geometry=tetrahedron("Si", "Cl", 2.0191),
     refs=[("R. R. Ryan, K. Hedberg, Effect of temperature on the structure of gaseous molecules. II. An " *
            "electron-diffraction investigation of the molecular structures of B2Cl4 and SiCl4, J. Chem. Phys. 50, " *
            "4986 (1969)", "10.1063/1.1670995")]),
    (name="Silicon tetrafluoride", formula="SiF4", category="Inorganic: nonmetal halides", pg="Td",
     geometry=tetrahedron("Si", "F", 1.554),
     refs=[("M. Takami, H. Kuze, Infrared-microwave double resonance spectroscopy of the SiF4 ν3 fundamental " *
            "using a tunable diode laser, J. Chem. Phys. 78, 2204 (1983)", "10.1063/1.445063")]),
    (name="Nitrogen trifluoride", formula="NF3", category="Inorganic: nonmetal halides", use="etch gas", pg="C3v",
     geometry=pyramid("N", "F", 1.3648, 102.37),
     refs=[("S. E. Novick, W. Chen, M. R. Munrow, K. J. Grant, Hyperfine structure in the microwave spectrum of " *
            "NF3, J. Mol. Spectrosc. 179, 219 (1996)", "10.1006/jmsp.1996.0199"), KUCHITSU]),
    (name="Carbon tetrafluoride", formula="CF4", category="Organic: halocarbons", use="etch gas", pg="Td",
     geometry=tetrahedron("C", "F", 1.3151), refs=[KUCHITSU]),
    # Etch and dopant gases
    (name="Boron trifluoride", formula="BF3", category="Inorganic: nonmetal halides", use="dopant source", pg="D3h",
     geometry=planar_xy3("B", "F", 1.30731),
     refs=[("A. G. Maki, T. Masiello, T. A. Blake, The high-resolution infrared spectrum of 10BF3 from 400 to " *
            "4600 cm-1, J. Mol. Struct. 742, 3 (2005)", "10.1016/j.molstruc.2004.11.089")]),
    (name="Boron trichloride", formula="BCl3", category="Inorganic: nonmetal halides", use="etch gas", pg="D3h",
     geometry=planar_xy3("B", "Cl", 1.7421),
     refs=[("S. Konaka, Y. Murata, K. Kuchitsu, Y. Morino, The molecular structure and force constants of boron " *
            "trifluoride and boron trichloride, Bull. Chem. Soc. Jpn. 39, 1134 (1966)", "10.1246/bcsj.39.1134")],
     note="Electron-diffraction r_g (thermal average), not an equilibrium r_e."),
    (name="Boron tribromide", formula="BBr3", category="Inorganic: nonmetal halides", pg="D3h",
     geometry=planar_xy3("B", "Br", 1.896),
     refs=[("K.-G. Martinsen, N. Vogt, H. V. Volden, V. S. Lyutsarev, J. Vogt, Molecular structure and force field " *
            "of boron tribromide as determined from combined analysis of gas electron diffraction and spectroscopic " *
            "data and supported by quantum-chemical density-functional calculations, J. Mol. Struct. 385, 159 (1996)",
            "10.1016/S0022-2860(96)09435-5")],
     note="Equilibrium estimate r_he from electron diffraction and spectroscopy (r_g = 1.900 Å)."),
    (name="Sulfur hexafluoride", formula="SF6", category="Inorganic: nonmetal halides", use="etch gas", pg="Oh",
     geometry=octahedron("S", "F", 1.5560),
     refs=[("V. Boudon, J. L. Doménech, D. Bermejo, H. Willner, High-resolution Raman spectroscopy of the ν1 region " *
            "and Raman-Raman double resonance spectroscopy of the 2ν1-ν1 band of 32SF6 and 34SF6. Determination of " *
            "the equilibrium bond length of sulfur hexafluoride, J. Mol. Spectrosc. 228, 392 (2004)",
            "10.1016/j.jms.2004.04.006")]),
    (name="Fluoromethane", formula="CH3F", category="Organic: halocarbons", use="etch gas", pg="C3v",
     geometry=xy3z("C", "H", "F", 1.087, 110.2, 1.383),
     refs=[("J. Demaison, J. Breidung, W. Thiel, D. Papoušek, The equilibrium structure of methyl fluoride, " *
            "Struct. Chem. 10, 129 (1999)", "10.1023/A:1022085314343")],
     note="Semi-experimental r_e: experimental rotational constants combined with CCSD(T) values by least squares."),
    (name="Hexafluoroethane", formula="C2F6", category="Organic: halocarbons", use="etch gas", pg="D3d",
     geometry=staggered("C", "F", 1.56, 1.32, 109.5),
     refs=[("D. A. Swick, I. L. Karle, Structure and internal motion of C2F6 and Si2Cl6 vapors, J. Chem. Phys. 23, " *
            "1499 (1955)", "10.1063/1.1742337")],
     note="Early electron-diffraction structure of low precision (C–C ± 0.03 Å, C–C–F ± 1.5°): relax before use."),
    (name="Chlorine trifluoride", formula="ClF3", category="Inorganic: nonmetal halides", use="etch gas", pg="C2v",
     geometry=zmatrix([("Cl",), ("F", 1, 1.598), ("F", 1, 1.698, 2, 87.483), ("F", 1, 1.698, 2, 87.483, 3, 180.0)]),
     refs=[("D. F. Smith, The microwave spectrum and structure of chlorine trifluoride, J. Chem. Phys. 21, 609 (1953)",
            "10.1063/1.1698976")],
     note="Planar T shape; effective ground-state structure from the moments of inertia."),
    (name="Xenon difluoride", formula="XeF2", category="Inorganic: nonmetal halides", use="etch gas", pg="D∞h",
     geometry=linear(("F", "Xe", "F"), (1.9791, 1.9791)),
     refs=[("N. J. Brassington, H. G. M. Edwards, D. A. Long, Pure rotational Raman spectrum of XeF2, J. Chem. " *
            "Soc., Faraday Trans. 2 74, 1208 (1978)", "10.1039/F29787401208")],
     note="Ground-state r_0, not an equilibrium r_e."),
    (name="Germane", formula="GeH4", category="Inorganic: hydrides", use="deposition precursor", pg="Td",
     geometry=tetrahedron("Ge", "H", 1.51714),
     refs=[("C. Richard, V. Boudon, A. Rizopoulos, J. Vander Auwera, F. Kwabia Tchana, Line positions and intensities " *
            "for the ν2/ν4 bands of 5 isotopologues of germane near 11.5 µm, J. Quant. Spectrosc. Radiat. Transfer " *
            "260, 107474 (2021)", "10.1016/j.jqsrt.2020.107474")]),
    (name="Phosphorus trifluoride", formula="PF3", category="Inorganic: nonmetal halides", pg="C3v",
     geometry=pyramid("P", "F", 1.56099, 97.567),
     refs=[("H. Najib, Experimental rovibrational constants and equilibrium structure of phosphorus trifluoride, " *
            "J. Mol. Spectrosc. 305, 17 (2014)", "10.1016/j.jms.2014.09.008")]),
    (name="Phosphorus pentafluoride", formula="PF5", category="Inorganic: nonmetal halides", pg="D3h",
     geometry=bipyramid("P", "F", 1.580, 1.532),
     refs=[("H. Kurimura, S. Yamamoto, T. Egawa, K. Kuchitsu, Geometrical structure of PF5 as studied by gas " *
            "electron diffraction and spectroscopic data, J. Mol. Struct. 140, 79 (1986)", "10.1016/0022-2860(86)80149-1")],
     note="Electron-diffraction r_g; the authors estimate r_e = 1.576 Å (axial) and 1.529 Å (equatorial)."),
    (name="Arsenic trifluoride", formula="AsF3", category="Inorganic: nonmetal halides", pg="C3v",
     geometry=pyramid("As", "F", 1.7041, 95.767),
     refs=[("J. G. Smith, The anharmonic force field of arsenic trifluoride from the millimetre wave spectra of " *
            "vibrationally excited states, Mol. Phys. 35, 461 (1978)", "10.1080/00268977800100341")]),
    # ALD/CVD precursors: metal alkyls and carbonyls
    (name="Trimethylborane", formula="B(CH3)3", category="Organometallic: metal alkyls", use="deposition precursor", pg="C3h",
     geometry=trimethyl("B", 1.5783, 1.1138, 111.9),
     refs=[("L. S. Bartell, B. L. Carroll, Electron-diffraction study of the structure of B(CH3)3, J. Chem. Phys. 42, " *
            "3076 (1965)", "10.1063/1.1696383")],
     note="Electron-diffraction mean distances. The methyl groups rotate freely; one orientation is chosen here " *
          "(one C–H of each methyl in the BC3 plane) with C–B–C = 120°."),
    (name="Trimethylaluminium", formula="Al(CH3)3", category="Organometallic: metal alkyls", use="deposition precursor", pg="C3h",
     geometry=trimethyl("Al", 1.957, 1.113, 111.7),
     refs=[("A. Almenningen, S. Halvorsen, A. Haaland, A gas phase electron diffraction investigation of the " *
            "molecular structures of trimethylaluminium monomer and dimer, Acta Chem. Scand. 25, 1937 (1971)",
            "10.3891/acta.chem.scand.25-1937")],
     note="The monomer (r_a). The methyl groups rotate freely; one orientation is chosen here (one C–H of each " *
          "methyl in the AlC3 plane). At room temperature the vapour is mostly the dimer, Al2(CH3)6."),
    (name="Trimethylaluminium dimer", formula="Al2(CH3)6", category="Organometallic: metal alkyls", use="deposition precursor", pg="C2v",
     geometry=bridged_trimethyl_dimer("Al", 2.619, 1.957, 2.140, 117.3, 1.117, 112.1),
     refs=[("A. Almenningen, S. Halvorsen, A. Haaland, A gas phase electron diffraction investigation of the " *
            "molecular structures of trimethylaluminium monomer and dimer, Acta Chem. Scand. 25, 1937 (1971)",
            "10.3891/acta.chem.scand.25-1937")],
     note="Electron diffraction (r_a): D2h Al2C6 skeleton with a mean C–H. The hydrogens follow the authors' model " *
          "(terminal methyls staggered); the bridging methyls rotate freely and one orientation is chosen here."),
    (name="Dimethylzinc", formula="Zn(CH3)2", category="Organometallic: metal alkyls", use="deposition precursor", pg="D3d",
     geometry=dimethyl("Zn", 1.930, 1.100, 112.5),
     refs=[("A. Almenningen, T. U. Helgaker, A. Haaland, S. Samdal, The molecular structures of dimethyl-, diethyl- " *
            "and dipropylzinc determined by gas phase electron diffraction. Normal coordinate analysis and ab initio " *
            "molecular orbital calculations on dimethylzinc, Acta Chem. Scand. A 36, 159 (1982)",
            "10.3891/acta.chem.scand.36a-0159")],
     note="Electron diffraction (r_a). The methyl groups rotate freely; the staggered form is chosen here."),
    (name="Trimethylindium", formula="In(CH3)3", category="Organometallic: metal alkyls", use="deposition precursor", pg="C3h",
     geometry=trimethyl("In", 2.161, 1.124, 109.8),
     refs=[("T. Fjeldberg, A. Haaland, R. Seip, Q. Shen, J. Weidlein, The molecular structures of trimethylindium " *
            "and trimethylthallium determined by gas electron diffraction, Acta Chem. Scand. A 36, 495 (1982)",
            "10.3891/acta.chem.scand.36a-0495")],
     note="Electron diffraction (r_a), the authors' model: planar InC3 and one C–H of each methyl in that plane."),
    (name="Nickel tetracarbonyl", formula="Ni(CO)4", category="Organometallic: metal carbonyls", use="deposition precursor", pg="Td",
     geometry=carbonyl("Ni", 4, 1.838, 1.141),
     refs=[("L. Hedberg, T. Iijima, K. Hedberg, Nickel tetracarbonyl, Ni(CO)4. I. Molecular structure by gaseous " *
            "electron diffraction. II. Refinement of quadratic force field, J. Chem. Phys. 70, 3224 (1979)",
            "10.1063/1.437911")],
     note="Electron-diffraction r_g. Highly toxic."),
    (name="Molybdenum hexacarbonyl", formula="Mo(CO)6", category="Organometallic: metal carbonyls", use="deposition precursor", pg="Oh",
     geometry=carbonyl("Mo", 6, 2.0632, 1.1450),
     refs=[("S. P. Arnesen, H. M. Seip, Studies on the failure of the first Born approximation in electron " *
            "diffraction. V. Molybdenum- and tungsten hexacarbonyl, Acta Chem. Scand. 20, 2711 (1966)",
            "10.3891/acta.chem.scand.20-2711")],
     note="Electron-diffraction r_g."),
    (name="Tungsten hexacarbonyl", formula="W(CO)6", category="Organometallic: metal carbonyls", use="deposition precursor", pg="Oh",
     geometry=carbonyl("W", 6, 2.0586, 1.1480),
     refs=[("S. P. Arnesen, H. M. Seip, Studies on the failure of the first Born approximation in electron " *
            "diffraction. V. Molybdenum- and tungsten hexacarbonyl, Acta Chem. Scand. 20, 2711 (1966)",
            "10.3891/acta.chem.scand.20-2711")],
     note="Electron-diffraction r_g."),
    # ALD/CVD precursors: halides
    (name="Titanium tetrachloride", formula="TiCl4", category="Inorganic: metal halides", use="deposition precursor", pg="Td",
     geometry=tetrahedron("Ti", "Cl", 2.185),
     refs=[("M. Kimura, K. Kimura, M. Aoki, S. Shibata, Investigation on the molecular structures of titanium " *
            "tetrachloride and zirconium tetrachloride by gas electron diffraction, Bull. Chem. Soc. Jpn. 29, 95 (1956)",
            "10.1246/bcsj.29.95")],
     note="Early electron-diffraction value (± 0.01 Å, distance type not stated)."),
    (name="Zirconium tetrachloride", formula="ZrCl4", category="Inorganic: metal halides", use="deposition precursor", pg="Td",
     geometry=tetrahedron("Zr", "Cl", 2.32),
     refs=[("M. Kimura, K. Kimura, M. Aoki, S. Shibata, Investigation on the molecular structures of titanium " *
            "tetrachloride and zirconium tetrachloride by gas electron diffraction, Bull. Chem. Soc. Jpn. 29, 95 (1956)",
            "10.1246/bcsj.29.95")],
     note="Early electron-diffraction value (± 0.01 Å, distance type not stated)."),
    (name="Tin tetrachloride", formula="SnCl4", category="Inorganic: metal halides", use="deposition precursor", pg="Td",
     geometry=tetrahedron("Sn", "Cl", 2.2808),
     refs=[("H. Fujii, M. Kimura, The molecular structure of tin tetrachloride as determined by gas-electron " *
            "diffraction, Bull. Chem. Soc. Jpn. 43, 1933 (1970)", "10.1246/bcsj.43.1933")],
     note="Electron-diffraction r_g."),
    (name="Molybdenum hexafluoride", formula="MoF6", category="Inorganic: metal halides", use="deposition precursor", pg="Oh",
     geometry=octahedron("Mo", "F", 1.820), refs=[("H. M. Seip, R. Seip, Studies on the failure of the first Born approximation in electron diffraction. " *
            "IV. Molybdenum- and tungsten hexafluoride, Acta Chem. Scand. 20, 2698 (1966)",
            "10.3891/acta.chem.scand.20-2698")], note="Electron-diffraction r_g."),
    (name="Tungsten hexafluoride", formula="WF6", category="Inorganic: metal halides", use="deposition precursor", pg="Oh",
     geometry=octahedron("W", "F", 1.832), refs=[("H. M. Seip, R. Seip, Studies on the failure of the first Born approximation in electron diffraction. " *
            "IV. Molybdenum- and tungsten hexafluoride, Acta Chem. Scand. 20, 2698 (1966)",
            "10.3891/acta.chem.scand.20-2698")], note="Electron-diffraction r_g."),
    (name="Tungsten hexachloride", formula="WCl6", category="Inorganic: metal halides", use="deposition precursor", pg="Oh",
     geometry=octahedron("W", "Cl", 2.281),
     refs=[("A. Haaland, K.-G. Martinsen, S. Shlykov, The molecular structure of tungsten hexachloride by gas " *
            "electron diffraction, Acta Chem. Scand. 46, 1208 (1992)", "10.3891/acta.chem.scand.46-1208")],
     note="Electron-diffraction r_a."),
    (name="Aluminium chloride dimer", formula="Al2Cl6", category="Inorganic: metal halides", use="deposition precursor", pg="D2h",
     geometry=bridged_halide_dimer("Al", "Cl", 2.250, 90.0, 2.061, 122.1),
     refs=[("K. Aarset, Q. Shen, H. Thomassen, A. D. Richardson, K. Hedberg, Molecular structure of the aluminum " *
            "halides, Al2Cl6, AlCl3, Al2Br6, AlBr3, and AlI3, obtained by gas-phase electron-diffraction and ab " *
            "initio molecular orbital calculations, J. Phys. Chem. A 103, 1644 (1999)", "10.1021/jp9842042")],
     note="Electron-diffraction r_g, planar Al2Cl2 ring (the thermal average is puckered); Cl(t)–Al–Cl(t) is " *
          "known only to ± 3°."),
    (name="Trichlorosilane", formula="SiHCl3", category="Inorganic: nonmetal halides", use="deposition precursor", pg="C3v",
     geometry=xy3z("Si", "Cl", "H", 2.0118, 110.60, 1.4655),
     refs=[("M. Mitzlaff, R. Holm, H. Hartmann, Mikrowellenspektrum, Struktur und Dipolmoment von Trichlorsilan " *
            "und Methyltrichlorsilan, Z. Naturforsch. A 22, 1415 (1967)", "10.1515/zna-1967-0922")],
     note="Substitution structure r_s from microwave spectra."),
    (name="Hexachlorodisilane", formula="Si2Cl6", category="Inorganic: nonmetal halides", use="deposition precursor", pg="D3d",
     geometry=staggered("Si", "Cl", 2.24, 2.01, 109.5),
     refs=[("D. A. Swick, I. L. Karle, Structure and internal motion of C2F6 and Si2Cl6 vapors, J. Chem. Phys. 23, " *
            "1499 (1955)", "10.1063/1.1742337")],
     note="Early electron-diffraction structure of low precision (Si–Si ± 0.06 Å): relax before use."),
    # Organic solvents and oxygen, nitrogen and sulfur compounds
    (name="Dimethyl sulfoxide", formula="(CH3)2SO", category="Organic: nitrogen and sulfur compounds", use="solvent",
     pg="Cs",
     geometry=cartesian([("S", 0.0, 0.1432, 0.4202), ("O", 0.0, 1.4024, -0.3667),
                         ("C", 1.3425, -0.8664, -0.2227), ("C", -1.3425, -0.8664, -0.2227),
                         ("H", 2.2553, -0.4311, 0.0751), ("H", 1.2255, -0.9383, -1.3113), ("H", 1.3057, -1.8579, 0.2349),
                         ("H", -2.2553, -0.4311, 0.0751), ("H", -1.2255, -0.9383, -1.3113),
                         ("H", -1.3057, -1.8579, 0.2349)]),
     refs=[("W. Feder, H. Dreizler, H. D. Rudolph, V. Typke, r_s-Struktur von Dimethylsulfoxid im Vergleich zur " *
            "r_0-Struktur, Z. Naturforsch. A 24, 266 (1969)", "10.1515/zna-1969-0215")],
     note="Substitution (r_s) coordinates of Table 3. The authors note that one hydrogen per methyl group (H I, " *
          "C–H 1.054 Å) is poorly located: relax before use."),
    (name="Ethylene oxide (oxirane)", formula="C2H4O", category="Organic: alcohols and ethers", pg="C2v",
     geometry=ring_ch2("O", 1.466, 1.431, 1.085, 116.61, 21.99),
     refs=[("C. Hirose, The microwave spectra and r_0, r_s, and r_m structures of ethylene oxide, Bull. Chem. Soc. " *
            "Jpn. 47, 1311 (1974)", "10.1246/bcsj.47.1311")],
     note="Substitution structure r_s."),
    (name="Ketene", formula="H2CCO", category="Organic: carbonyl compounds", pg="C2v",
     geometry=([("C"), "C", "O", "H", "H"],
               [[0.0, 0, 0], [1.31212, 0, 0], [1.31212 + 1.16030, 0, 0],
                [-1.07576 * cosd(121.781 / 2), 1.07576 * sind(121.781 / 2), 0],
                [-1.07576 * cosd(121.781 / 2), -1.07576 * sind(121.781 / 2), 0]]),
     refs=[("A. L. L. East, W. D. Allen, S. J. Klippenstein, The anharmonic force field and equilibrium molecular " *
            "structure of ketene, J. Chem. Phys. 102, 8506 (1995)", "10.1063/1.468842")],
     note="Semi-experimental r_e: rotational constants of six isotopologues with ab initio vibrational corrections."),
    (name="Formamide", formula="HCONH2", category="Organic: carbonyl compounds", pg="C1",
     geometry=zmatrix([("C",), ("O", 1, 1.193), ("N", 1, 1.376, 2, 123.8), ("H", 1, 1.102, 3, 113.233, 2, 180.0),
                       ("H", 3, 1.014, 1, 117.15, 2, FORMAMIDE_TWIST[1]),
                       ("H", 3, 1.002, 1, 120.617, 4, FORMAMIDE_TWIST[2])]),
     refs=[("C. C. Costain, J. M. Dowling, Microwave spectrum and molecular structure of formamide, J. Chem. Phys. " *
            "32, 158 (1960)", "10.1063/1.1700891")],
     note="The authors' slightly pyramidal NH2 (dihedrals 7° and 12°, ± 5°); later work finds formamide " *
          "effectively planar."),
    (name="Methyl formate", formula="HCOOCH3", category="Organic: carbonyl compounds", pg="Cs",
     geometry=methyl_formate(),
     refs=[("R. F. Curl, Microwave spectrum, barrier to internal rotation, and structure of methyl formate, " *
            "J. Chem. Phys. 30, 1529 (1959)", "10.1063/1.1730232")],
     note="cis conformer (methyl cis to C=O). The author assumed the methyl axis along the O–CH3 bond; the methyl " *
          "hydrogens stagger the C–O bond (one in the heavy-atom plane)."),
    (name="Chloromethane", formula="CH3Cl", category="Organic: halocarbons", pg="C3v",
     geometry=xy3z_axial("C", "H", "Cl", 1.083484, 108.3723, 1.777716),
     refs=[("M. Melosso, P. Stoppa, D. Alvarado-Jiménez, F. Tamassia, C. Sapienza, L. Bizzocchi, L. Dore, " *
            "C. Puzzarini, A. Pietropolli Charmet, N. Tasinato, Completing the spectral mosaic of chloromethane by " *
            "adding the CHD2Cl missing piece through the interplay of rotational/vibrational spectroscopy and " *
            "quantum chemical calculations, Molecules 30, 1604 (2025)", "10.3390/molecules30071604")],
     note="Semi-experimental r_e (Table 5): rotational constants with CCSD(T) vibrational corrections."),
    (name="Dichloromethane", formula="CH2Cl2", category="Organic: halocarbons", use="solvent", pg="C2v",
     geometry=xy2z2("C", "H", "Cl", 1.082, 112.967, 1.7724, 111.783),
     refs=[("R. J. Myers, W. D. Gwinn, The microwave spectra, structure, dipole moment, and chlorine nuclear " *
            "quadrupole coupling constants of methylene chloride, J. Chem. Phys. 20, 1420 (1952)", "10.1063/1.1700773")],
     note="Effective (r_0) structure; C–H and H–C–H are the authors' values extrapolated to infinite hydrogen mass."),
    (name="Chloroform", formula="CHCl3", category="Organic: halocarbons", use="solvent", pg="C3v",
     geometry=xy3z_axial("C", "Cl", "H", 1.760, 108.23, 1.080),
     refs=[("J.-M. Colmont, D. Priem, P. Dréan, J. Demaison, J. E. Boggs, Rotational spectra of the isotopic species " *
            "of chloroform: experimental and ab initio structures, J. Mol. Spectrosc. 191, 158 (1998)",
            "10.1006/jmsp.1998.7623")],
     note="Approximate r_e combining the experiment with ab initio calculations."),
    (name="Carbon tetrachloride", formula="CCl4", category="Organic: halocarbons", pg="Td",
     geometry=tetrahedron("C", "Cl", 1.760),
     refs=[("L. Pauling, L. O. Brockway, A study of the methods of interpretation of electron-diffraction " *
            "photographs of gas molecules, with results for benzene and carbon tetrachloride, J. Chem. Phys. 2, 867 " *
            "(1934)", "10.1063/1.1749410")],
     note="An early (1934) electron-diffraction value, ± 0.005 Å: relax before use."),
    (name="Acetonitrile", formula="CH3CN", category="Organic: nitrogen and sulfur compounds", use="solvent", pg="C3v",
     geometry=methyl_nitrile(1.468, 1.159, 1.107, 109.7),
     refs=[("K. Karakida, T. Fukuyama, K. Kuchitsu, Molecular structures of hydrogen cyanide and acetonitrile as " *
            "studied by gas electron diffraction, Bull. Chem. Soc. Jpn. 47, 299 (1974)", "10.1246/bcsj.47.299")],
     note="Electron-diffraction r_g bond lengths (C–H longer than r_e) and r_av angle."),
    (name="Pyridine", formula="C5H5N", category="Organic: heterocycles", use="solvent", pg="C2v",
     geometry=pyridine(1.3355, 1.3905, 1.3883, 116.973, 123.748, 118.544, (1.0821, 1.0801, 1.0806),
                       (120.316, 121.332, 120.778)),
     refs=[("M. A. Zdanovskaia, B. J. Esselman, S. M. Kougias, M. G. Atwood, G. H. Jones, J. F. Stanton, R. C. Woods, " *
            "R. J. McMahon, Precise semi-experimental equilibrium (r_e^SE) structure of pyridine from 32 " *
            "isotopologues, J. Phys. Chem. A 129, 10737 (2025)", "10.1021/acs.jpca.5c07184")],
     note="Semi-experimental r_e (Table 4)."),
    (name="Furan", formula="C4H4O", category="Organic: heterocycles", pg="C2v",
     geometry=ring5("O", 1.3595, 1.3547, 106.576, 110.720, 1.0736, 115.868, 1.0749, 126.352),
     refs=[("K. T. Mehalko, B. J. Esselman, P. R. Franke, M. S. Holdren, G. Wenzel, D. A. Stewart, Z. P. Burke, " *
            "B. A. McGuire, R. C. Woods, R. J. McMahon, Structure of furan from 44 isotopologues: resolving " *
            "differences in the theoretical (r_e) and semi-experimental equilibrium (r_e^SE) structures, J. Phys. " *
            "Chem. A 130, 7783 (2026)", "10.1021/acs.jpca.6c04036")],
     note="Semi-experimental r_e, the authors' recommended values (Table 4)."),
    (name="Thiophene", formula="C4H4S", category="Organic: heterocycles", pg="C2v",
     geometry=ring5("S", 1.71049, 1.36564, 92.047, 111.608, 1.07714, 120.065, 1.07856, 123.414),
     refs=[("V. L. Orr, Y. Ichikawa, A. R. Patel, S. M. Kougias, K. Kobayashi, J. F. Stanton, B. J. Esselman, " *
            "R. C. Woods, R. J. McMahon, Precise equilibrium structure determination of thiophene (c-C4H4S) by " *
            "rotational spectroscopy: structure of a five-membered heterocycle containing a third-row atom, " *
            "J. Chem. Phys. 154, 244310 (2021)", "10.1063/5.0055267")],
     note="Semi-experimental r_e (Table IV)."),
    (name="Propyne", formula="CH3CCH", category="Organic: hydrocarbons", pg="C3v",
     geometry=propyne(1.45884, 1.20460, 1.08852, 110.591, 1.06146),
     refs=[("H. S. P. Müller, S. Thorwirth, F. Lewen, Rotational spectroscopy of singly 13C substituted isotopomers " *
            "of propyne and determination of a semi-empirical equilibrium structure, J. Mol. Struct. 1207, 127769 " *
            "(2020)", "10.1016/j.molstruc.2020.127769")],
     note="Semi-experimental r_e (Table 3, ae-CCSD(T) corrections), read in the authors' preprint arXiv:2001.11907."),
    (name="Allene", formula="H2CCCH2", category="Organic: hydrocarbons", pg="D2d",
     geometry=allene(1.312, 1.083, 120.6),
     refs=[("A. Almenningen, O. Bastiansen, M. Trætteberg, An electron diffraction investigation of the molecular " *
            "structure of allene, Acta Chem. Scand. 13, 1699 (1959)", "10.3891/acta.chem.scand.13-1699")],
     note="Early electron diffraction without stated uncertainties: C–H 1.083 Å and H–C=C 120.6° from the fit, " *
          "C=C 1.312 Å as the authors' best value (1.308 Å from the fit)."),
    (name="Cyclopropane", formula="C3H6", category="Organic: hydrocarbons", pg="D3h",
     geometry=cyclopropane(1.5030, 1.0786, 114.97),
     refs=[("J. Gauss, D. Cremer, J. F. Stanton, The r_e structure of cyclopropane, J. Phys. Chem. A 104, 1319 (2000)",
            "10.1021/jp993189d")],
     note="Semi-experimental r_e: rotational constants with computed vibrational corrections."),
    (name="Cyclobutane", formula="C4H8", category="Organic: hydrocarbons", pg="D2d",
     geometry=cyclobutane(1.552, 1.093, 106.4, 27.9, 6.2),
     refs=[("T. Egawa, T. Fukuyama, S. Yamamoto, F. Takabayashi, H. Kambara, T. Ueda, K. Kuchitsu, Molecular " *
            "structure and puckering potential function of cyclobutane studied by gas electron diffraction and " *
            "infrared spectroscopy, J. Chem. Phys. 86, 6018 (1987)", "10.1063/1.452489")],
     note="r_z structure (electron diffraction with infrared data): puckered ring, axial and equatorial C–H equal."),
    # Radicals and fragments (surface science, plasma and CVD intermediates)
    (name="Methylene", formula="CH2", category="Organic: radicals", pg="C2v",
     geometry=zmatrix([("C",), ("H", 1, 1.0748), ("H", 1, 1.0748, 2, 133.84)]),
     refs=[("P. R. Bunker, P. Jensen, A refined potential surface for the X̃ 3B1 electronic state of methylene CH2, " *
            "J. Chem. Phys. 79, 1224 (1983)", "10.1063/1.445927")],
     note="Equilibrium structure of the X̃ 3B1 ground state, fitted to experimental energy levels. " *
          spin_note("triplet")),
    (name="Formyl radical", formula="HCO", category="Organic: radicals", pg="Cs",
     geometry=zmatrix([("C",), ("O", 1, 1.17), ("H", 1, 1.11, 2, 127.0)]),
     refs=[("J. A. Austin, D. H. Levy, C. A. Gottlieb, H. E. Radford, Microwave spectrum of the HCO radical, " *
            "J. Chem. Phys. 60, 207 (1974)", "10.1063/1.1680769")],
     note=odd_electron_note("Low-precision geometry (no uncertainties given): relax before use.")),
    (name="Trifluoromethyl radical", formula="CF3", category="Organic: radicals", pg="C3v",
     geometry=pyramid("C", "F", 1.318, 110.76),
     refs=[("C. Yamada, E. Hirota, Infrared diode laser spectroscopy of the CF3 ν3 band, J. Chem. Phys. 78, 1703 " *
            "(1983)", "10.1063/1.444969")],
     note="Pyramidal r_0 structure. " * odd_electron_note()),
    (name="Amidogen radical", formula="NH2", category="Inorganic: radicals", pg="C2v",
     geometry=zmatrix([("N",), ("H", 1, 1.0245), ("H", 1, 1.0245, 2, 103.34)]),
     refs=[("I. Morino, K. Kawaguchi, Fourier transform far-infrared spectroscopy of the NH2, NHD, and ND2 " *
            "radicals, J. Mol. Spectrosc. 182, 428 (1997)", "10.1006/jmsp.1996.7221")],
     note="r_0 structure. " * odd_electron_note()),
    (name="Silylene", formula="SiH2", category="Inorganic: radicals", pg="C2v",
     geometry=zmatrix([("Si",), ("H", 1, 1.516), ("H", 1, 1.516, 2, 92.083)]),
     refs=[("I. Dubois, The absorption spectrum of the free SiH2 radical, Can. J. Phys. 46, 2485 (1968)",
            "10.1139/p68-608")],
     note="Ground state 1A1 (singlet, closed shell) from the electronic spectrum; no uncertainties given."),
    (name="Silyl radical", formula="SiH3", category="Inorganic: radicals", pg="C3v",
     geometry=pyramid("Si", "H", 1.468, 110.5),
     refs=[("C. Yamada, E. Hirota, Detection of the silyl radical SiH3 by infrared diode-laser spectroscopy, Phys. " *
            "Rev. Lett. 56, 923 (1986)", "10.1103/PhysRevLett.56.923")],
     note="Pyramidal effective structure (no uncertainties given). " * odd_electron_note()),
    (name="Hydroperoxyl radical", formula="HO2", category="Inorganic: radicals", pg="Cs",
     geometry=zmatrix([("O",), ("O", 1, 1.33054), ("H", 1, 0.9707, 2, 104.29)]),
     refs=[("K. G. Lubic, T. Amano, H. Uehara, K. Kawaguchi, E. Hirota, The ν1 band of the DO2 radical by difference " *
            "frequency laser and diode laser spectroscopy: the equilibrium structure of the hydroperoxyl radical, " *
            "J. Chem. Phys. 81, 4826 (1984)", "10.1063/1.447508")],
     note="Equilibrium structure r_e. " * odd_electron_note()),
    (name="Nitrogen dioxide", formula="NO2", category="Inorganic: oxides and sulfides", pg="C2v",
     geometry=zmatrix([("N",), ("O", 1, 1.1949), ("O", 1, 1.1949, 2, 133.75)]),
     refs=[("J. Demaison, J. Liévin, A. Perrin, Equilibrium structure of nitrogen dioxide revisited, Mol. Phys. 122, " *
            "e2318016 (2024)", "10.1080/00268976.2024.2318016")],
     note="Experimental r_e (Table 1) of the X̃ 2A1 ground state. " * spin_note("doublet")),
    (name="Methylidyne radical", formula="CH", category="Organic: radicals", pg="C∞v",
     geometry=linear(("C", "H"), (1.11983,)),
     refs=[("P. F. Bernath, The vibration-rotation emission spectrum of CH(X 2Π), J. Chem. Phys. 86, 4838 (1987)",
            "10.1063/1.452677")],
     note="Equilibrium r_e of the X 2Π ground state. " * spin_note("doublet")),
    (name="Mercapto radical", formula="SH", category="Inorganic: radicals", pg="C∞v",
     geometry=linear(("S", "H"), (1.3503,)),
     refs=[("D. A. Ramsay, Absorption spectra of SH and SD produced by flash photolysis of H2S and D2S, J. Chem. " *
            "Phys. 20, 1920 (1952)", "10.1063/1.1700339")],
     note="Ground-state r_0 (2Π), no uncertainty given. " * spin_note("doublet")),
    (name="Chlorine monoxide", formula="ClO", category="Inorganic: radicals", pg="C∞v",
     geometry=linear(("Cl", "O"), (1.569,)),
     refs=[("T. Amano, S. Saito, E. Hirota, Y. Morino, D. R. Johnson, F. X. Powell, The microwave spectrum of the ClO " *
            "radical, Ohio State University International Symposium on Molecular Spectroscopy, 1968, paper E2, " *
            "http://hdl.handle.net/1811/15452", "")],
     note="Estimated equilibrium distance (± 0.001 Å) from the authors' symposium abstract; 2Π ground state. " *
          spin_note("doublet")),
    (name="Silicon monoxide", formula="SiO", category="Inorganic: oxides and sulfides", pg="C∞v",
     geometry=linear(("Si", "O"), (1.50973,)),
     refs=[("T. Törring, Das Mikrowellenrotationsspektrum des SiO, Z. Naturforsch. A 23, 777 (1968)",
            "10.1515/zna-1968-0527")],
     note="Equilibrium r_e from the microwave spectrum (Table 2)."),
    # Oxoacids and more oxides
    (name="Nitric acid", formula="HNO3", category="Inorganic: oxoacids", use="wet etch and cleaning", pg="Cs",
     geometry=zmatrix([("N",), ("O", 1, 1.406), ("O", 1, 1.211, 2, 115.883), ("O", 1, 1.199, 2, 113.85, 3, 180.0),
                       ("H", 2, 0.964, 1, 102.15, 3, 0.0)]),
     refs=[("A. P. Cox, J. M. Riveros, Microwave spectrum and structure of nitric acid, J. Chem. Phys. 42, 3106 (1965)",
            "10.1063/1.1696387")],
     note="Planar substitution (r_s) structure; the O–H points towards the cis N–O."),
    (name="Nitrous acid (trans)", formula="HONO", category="Inorganic: oxoacids", pg="Cs",
     geometry=zmatrix([("N",), ("O", 1, 1.432), ("O", 1, 1.170, 2, 110.7), ("H", 2, 0.958, 1, 102.1, 3, 180.0)]),
     refs=[("A. P. Cox, A. H. Brittain, D. J. Finnigan, Microwave spectrum, structure, dipole moment and quadrupole " *
            "coupling constants of cis and trans nitrous acids, Trans. Faraday Soc. 67, 2179 (1971)",
            "10.1039/tf9716702179")],
     note="The planar trans isomer (H–O–N=O 180°); no uncertainties given."),
    (name="Hypofluorous acid", formula="HOF", category="Inorganic: oxoacids", pg="Cs",
     geometry=zmatrix([("O",), ("F", 1, 1.442), ("H", 1, 0.964, 2, 97.2)]),
     refs=[("H. Kim, E. F. Pearson, E. H. Appelman, Millimeter-wave spectrum and structure of hypofluorous acid: HOF " *
            "and DOF, J. Chem. Phys. 56, 1 (1972)", "10.1063/1.1676830")],
     note="Effective (r_0) structure from HOF and DOF."),
    (name="Sulfuric acid", formula="H2SO4", category="Inorganic: oxoacids", use="wet etch and cleaning", pg="C2",
     geometry=sulfuric_acid(1.574, 1.422, 0.97, 101.3, 123.3, 88.4, 108.5, -90.9),
     refs=[("R. L. Kuczkowski, R. D. Suenram, F. J. Lovas, Microwave spectrum, structure, and dipole moment of " *
            "sulfuric acid, J. Am. Chem. Soc. 103, 2561 (1981)", "10.1021/ja00400a013")],
     note="The C2 conformer from the microwave spectrum (structure type not stated)."),
    (name="Sulfur trioxide", formula="SO3", category="Inorganic: oxides and sulfides", pg="D3h",
     geometry=planar_xy3("S", "O", 1.4175),
     refs=[("V. Meyer, D. H. Sutter, H. Dreizler, The centrifugally induced pure rotational spectrum and the structure " *
            "of sulfur trioxide. A microwave Fourier transform study of a nonpolar molecule, Z. Naturforsch. A 46, " *
            "710 (1991)", "10.1515/zna-1991-0811")],
     note="Equilibrium r_e (no uncertainty given; r_0 = 1.4198 Å)."),
    (name="Dinitrogen tetroxide", formula="N2O4", category="Inorganic: oxides and sulfides", pg="D2h",
     geometry=zmatrix([("N",), ("N", 1, 1.782), ("O", 1, 1.190, 2, 112.3), ("O", 1, 1.190, 2, 112.3, 3, 180.0),
                       ("O", 2, 1.190, 1, 112.3, 3, 0.0), ("O", 2, 1.190, 1, 112.3, 3, 180.0)]),
     refs=[("B. W. McClelland, G. Gundersen, K. Hedberg, Reinvestigation of the structure of dinitrogen tetroxide, " *
            "N2O4, by gaseous electron diffraction, J. Chem. Phys. 56, 4541 (1972)", "10.1063/1.1677901")],
     note="Planar electron-diffraction r_a structure (O–N–O 135.4°). The N–N bond, 1.782 Å, is unusually long " *
          "and beyond the default bond limit, so it is not drawn; ViewSettings(bond_cutoffs=Dict((:N, :N) => 1.85)) " *
          "shows it."),
    (name="Dichlorine monoxide", formula="Cl2O", category="Inorganic: oxides and sulfides", pg="C2v",
     geometry=zmatrix([("O",), ("Cl", 1, 1.69587), ("Cl", 1, 1.69587, 2, 110.886)]),
     refs=[("M. Nakata, M. Sugie, H. Takeo, C. Matsumura, T. Fukuyama, K. Kuchitsu, Structure of dichlorine monoxide " *
            "as studied by microwave spectroscopy. Determination of equilibrium structure by a modified mass " *
            "dependence method, J. Mol. Spectrosc. 86, 241 (1981)", "10.1016/0022-2852(81)90122-3")],
     note="Equilibrium r_e estimated from r_m structures."),
    (name="Carbon disulfide", formula="CS2", category="Inorganic: oxides and sulfides", use="solvent", pg="D∞h",
     geometry=linear(("S", "C", "S"), (1.55448, 1.55448)),
     refs=[("D. S. Kummli, H. M. Frey, S. Leutwyler, Femtosecond degenerate four-wave mixing of carbon disulfide: " *
            "high-accuracy rotational constants, J. Chem. Phys. 124, 144307 (2006)", "10.1063/1.2186642")],
     note="Semi-experimental r_e (measured rotational constant with ab initio corrections)."),
    (name="Oxygen difluoride", formula="OF2", category="Inorganic: nonmetal halides", pg="C2v",
     geometry=zmatrix([("O",), ("F", 1, 1.3896), ("F", 1, 1.3896, 2, 104.163)]),
     refs=[("A. R. Hilton, A. W. Jache, J. B. Beal, W. D. Henderson, R. J. Robinson, Millimeter wave spectrum and " *
            "molecular structure of oxygen difluoride, J. Chem. Phys. 34, 1137 (1961)", "10.1063/1.1731711")],
     note="An early effective structure from one isotopologue, no uncertainties given."),
    # Oxyhalides and molecular elements
    (name="Phosgene", formula="COCl2", category="Inorganic: oxyhalides", pg="C2v",
     geometry=zmatrix([("C",), ("O", 1, 1.1766), ("Cl", 1, 1.7365, 2, 180 - 111.91 / 2),
                       ("Cl", 1, 1.7365, 2, 180 - 111.91 / 2, 3, 180.0)]),
     refs=[("S. Yamamoto, M. Nakata, K. Kuchitsu, Equilibrium structure and anharmonic potential constants of phosgene " *
            "derived from rotational constants and electron diffraction intensity, J. Mol. Spectrosc. 112, 173 (1985)",
            "10.1016/0022-2852(85)90202-4")],
     note="Equilibrium r_e. Highly toxic."),
    (name="Carbonyl fluoride", formula="COF2", category="Inorganic: oxyhalides", pg="C2v",
     geometry=zmatrix([("C",), ("O", 1, 1.1717), ("F", 1, 1.3157, 2, 180 - 107.71 / 2),
                       ("F", 1, 1.3157, 2, 180 - 107.71 / 2, 3, 180.0)]),
     refs=[("M. Nakata, K. Kohata, T. Fukuyama, K. Kuchitsu, C. J. Wilkins, Molecular structure of carbonyl fluoride " *
            "as studied by gas electron diffraction and microwave data, J. Mol. Struct. 68, 271 (1980)",
            "10.1016/0022-2860(80)80278-X")],
     note="r_z structure from electron diffraction with microwave data."),
    (name="Thionyl chloride", formula="SOCl2", category="Inorganic: oxyhalides", pg="Cs",
     geometry=xyz2_pyramid("S", "O", "Cl", 1.435, 2.072, 108.00, 97.15),
     refs=[("S. Suzuki, M. Yamaguchi, M. Onda, T. Sakaizumi, O. Ohashi, I. Yamaguchi, Microwave spectroscopic study on " *
            "the molecular structure and the quadrupole coupling constants of thionyl chloride, J. Mol. Struct. 73, 41 " *
            "(1981)", "10.1016/0022-2860(81)85044-2")],
     note="r_0 structure derived with some assumptions (not stated in the abstract)."),
    (name="Sulfuryl chloride", formula="SO2Cl2", category="Inorganic: oxyhalides", pg="C2v",
     geometry=xy2z2("S", "O", "Cl", 1.41347, 123.129, 2.01124, 100.126),
     refs=[("J. Chen, Rotational and vibrational spectroscopy of sulphuryl chloride, sulphuryl fluoride and dimethyl " *
            "sulphone, Ph.D. thesis, Memorial University of Newfoundland (1991), " *
            "https://hdl.handle.net/20.500.14783/6953", "")],
     note="Ground-state average structure of 35Cl2 32S 16O2 from the thesis abstract."),
    (name="Phosphoryl chloride", formula="POCl3", category="Inorganic: oxyhalides", use="dopant source", pg="C3v",
     geometry=xy3z("P", "Cl", "O", 1.99, 103.6, 1.45),
     refs=[("D. Williams, J. Sheridan, W. Gordy, Microwave spectra and molecular structures of POF3, PSF3, POCl3, and " *
            "PSCl3, J. Chem. Phys. 20, 164 (1952)", "10.1063/1.1700162")],
     note="Early structure of low precision (± 0.03 Å, ± 2°): relax before use."),
    (name="Nitrosyl fluoride", formula="NOF", category="Inorganic: oxyhalides", pg="Cs",
     geometry=zmatrix([("N",), ("O", 1, 1.136), ("F", 1, 1.512, 2, 110.083)]),
     refs=[("K. S. Buckton, A. C. Legon, D. J. Millen, Microwave spectrum and structure of nitrosyl fluoride, Trans. " *
            "Faraday Soc. 65, 1975 (1969)", "10.1039/tf9696501975")],
     note="From three isotopic species (structure type not stated)."),
    (name="White phosphorus", formula="P4", category="Inorganic: elements", pg="Td",
     geometry=tetrahedral_x4("P", 2.1994),
     refs=[("B. M. Cossairt, C. C. Cummins, A. R. Head, D. L. Lichtenberger, R. J. F. Berger, S. A. Hayes, " *
            "N. W. Mitzel, G. Wu, On the molecular and electronic structures of AsP3 and P4, J. Am. Chem. Soc. 132, " *
            "8459 (2010)", "10.1021/ja102580d")],
     note="Gas electron diffraction (r_g)."),
    (name="Cyclo-octasulfur", formula="S8", category="Inorganic: elements", pg="D4d",
     geometry=crown_ring("S", 8, 2.055, 108.2),
     refs=[("S. J. Rettig, J. Trotter, Refinement of the structure of orthorhombic sulfur, α-S8, Acta Crystallogr. C 43, " *
            "2260 (1987)", "10.1107/S0108270187088152")],
     note="Crown ring from the mean crystal values (S–S corrected for libration, S–S–S 108.2°); the torsion " *
          "follows, 98.5° in the crystal."),
    (name="Selenium dioxide", formula="SeO2", category="Inorganic: oxides and sulfides", pg="C2v",
     geometry=zmatrix([("Se",), ("O", 1, 1.6076), ("O", 1, 1.6076, 2, 113.833)]),
     refs=[("H. Takeo, E. Hirota, Y. Morino, Equilibrium structure and potential function of selenium dioxide by " *
            "microwave spectroscopy, J. Mol. Spectrosc. 34, 370 (1970)", "10.1016/0022-2852(70)90020-2")],
     note="Equilibrium r_e of the gas-phase monomer."),
    # Main-group halides
    (name="Aluminium trichloride (monomer)", formula="AlCl3", category="Inorganic: metal halides",
     use="deposition precursor", pg="D3h",
     geometry=planar_xy3("Al", "Cl", 2.062),
     refs=[("K. Aarset, Q. Shen, H. Thomassen, A. D. Richardson, K. Hedberg, Molecular structure of the aluminum " *
            "halides, Al2Cl6, AlCl3, Al2Br6, AlBr3, and AlI3, obtained by gas-phase electron-diffraction and ab " *
            "initio molecular orbital calculations, J. Phys. Chem. A 103, 1644 (1999)", "10.1021/jp9842042")],
     note="Electron-diffraction r_g of the monomer (D3h assumed); the vapour is mostly the dimer Al2Cl6 at lower " *
          "temperature."),
    (name="Gallium trichloride (monomer)", formula="GaCl3", category="Inorganic: metal halides",
     use="deposition precursor", pg="D3h", geometry=planar_xy3("Ga", "Cl", 2.108), refs=[("A. Haaland, A. Hammel, K.-G. Martinsen, J. Tremmel, H. V. Volden, Molecular structures of monomeric gallium " *
            "trichloride, indium trichloride and lead tetrachloride by gas electron diffraction, J. Chem. Soc., Dalton " *
            "Trans. 1992, 2209", "10.1039/DT9920002209")],
     note="Electron-diffraction r_a of the monomer (about 380 °C)."),
    (name="Indium trichloride (monomer)", formula="InCl3", category="Inorganic: metal halides",
     use="deposition precursor", pg="D3h", geometry=planar_xy3("In", "Cl", 2.289), refs=[("A. Haaland, A. Hammel, K.-G. Martinsen, J. Tremmel, H. V. Volden, Molecular structures of monomeric gallium " *
            "trichloride, indium trichloride and lead tetrachloride by gas electron diffraction, J. Chem. Soc., Dalton " *
            "Trans. 1992, 2209", "10.1039/DT9920002209")],
     note="Electron-diffraction r_a of the monomer (about 480 °C)."),
    (name="Beryllium dichloride", formula="BeCl2", category="Inorganic: metal halides", pg="D∞h",
     geometry=linear(("Cl", "Be", "Cl"), (1.791, 1.791)),
     refs=[("A. G. Girichev, N. I. Giricheva, N. Vogt, G. V. Girichev, J. Vogt, Structural investigation of molecules " *
            "in the vapour over beryllium dichloride using electron diffraction and mass spectrometric data, J. Mol. " *
            "Struct. 384, 175 (1996)", "10.1016/S0022-2860(96)09333-7")],
     note="Equilibrium r_e^anh of the linear monomer (r_g = 1.798 Å). Highly toxic."),
    (name="Phosphorus trichloride", formula="PCl3", category="Inorganic: nonmetal halides", pg="C3v",
     geometry=pyramid("P", "Cl", 2.043, 100.10), refs=[("P. Kisliuk, C. H. Townes, The microwave spectra and molecular structure of phosphorus and arsenic " *
            "trichloride, J. Chem. Phys. 18, 1109 (1950)", "10.1063/1.1747872")],
     note="Early microwave structure (effectively r_0)."),
    (name="Arsenic trichloride", formula="AsCl3", category="Inorganic: nonmetal halides", pg="C3v",
     geometry=pyramid("As", "Cl", 2.161, 98.42), refs=[("P. Kisliuk, C. H. Townes, The microwave spectra and molecular structure of phosphorus and arsenic " *
            "trichloride, J. Chem. Phys. 18, 1109 (1950)", "10.1063/1.1747872")],
     note="Early microwave structure (effectively r_0)."),
    (name="Antimony trichloride", formula="SbCl3", category="Inorganic: nonmetal halides", pg="C3v",
     geometry=pyramid("Sb", "Cl", 2.333, 97.2),
     refs=[("S. Konaka, M. Kimura, The molecular structure of antimony trichloride as determined by gas-electron " *
            "diffraction, Bull. Chem. Soc. Jpn. 46, 404 (1973)", "10.1246/bcsj.46.404")],
     note="Electron-diffraction r_g with the angle θ_α."),
    (name="Nitrogen trichloride", formula="NCl3", category="Inorganic: nonmetal halides", pg="C3v",
     geometry=pyramid("N", "Cl", 1.759, 107.1),
     refs=[("H. B. Bürgi, D. Stedman, L. S. Bartell, Molecular structure of nitrogen trichloride as determined by " *
            "electron diffraction, J. Mol. Struct. 10, 31 (1971)", "10.1016/0022-2860(71)87058-8")],
     note="Electron-diffraction r_g. Explosive."),
    (name="Phosphorus pentachloride", formula="PCl5", category="Inorganic: nonmetal halides", pg="D3h",
     geometry=bipyramid("P", "Cl", 2.124, 2.020),
     refs=[("W. J. Adams, L. S. Bartell, Structure and vibrational assignments for PCl5. An electron diffraction " *
            "study, J. Mol. Struct. 8, 23 (1971)", "10.1016/0022-2860(71)80038-8")],
     note="Electron diffraction (Table 2): axial bonds 0.104 Å longer than equatorial."),
    # Silicon compounds
    (name="Methylsilane", formula="CH3SiH3", category="Organometallic: metal alkyls", use="deposition precursor",
     pg="C3v", geometry=methyl_xh3("Si", 1.8668, 1.093, 1.485, 107.667, 108.25),
     refs=[("R. W. Kilb, L. Pierce, Microwave spectrum, structure, and internal barrier of methyl silane, J. Chem. " *
            "Phys. 27, 108 (1957)", "10.1063/1.1743646")],
     note="Ground-state (effectively r_0) structure; methyl and silyl groups staggered, as the authors established."),
    (name="Chlorosilane", formula="SiH3Cl", category="Inorganic: nonmetal halides", pg="C3v",
     geometry=xy3z_axial("Si", "H", "Cl", 1.469, 108.43, 2.046),
     refs=[("G. Cazzoli, C. Puzzarini, J. Gauss, Rotational spectrum of silyl chloride: hyperfine structure and " *
            "equilibrium geometry, Mol. Phys. 110, 2359 (2012)", "10.1080/00268976.2012.680518")],
     note="Recommended equilibrium structure (semi-experimental r_e in agreement with theory)."),
    (name="Fluorosilane", formula="SiH3F", category="Inorganic: nonmetal halides", pg="C3v",
     geometry=xy3z_axial("Si", "H", "F", 1.4698, 108.29, 1.5906),
     refs=[("C. Puzzarini, G. Cazzoli, J. Gauss, Silyl fluoride: Lamb-dip spectra and equilibrium structure, " *
            "64th International Symposium on Molecular Spectroscopy, Columbus (2009), talk, " *
            "http://hdl.handle.net/1811/38302", ""),
           ("C. Puzzarini, G. Cazzoli, J. Gauss, J. Mol. Spectrosc. 262, 37 (2010) (the paper of this work)",
            "10.1016/j.jms.2010.05.001")],
     note="Semi-experimental r_e, read in the authors' symposium slides (the journal paper was not readable)."),
    (name="Disiloxane", formula="(SiH3)2O", category="Inorganic: oxides and sulfides", pg="C2v",
     geometry=disiloxane("Si", 1.634, 1.486, 144.1, 109.9),
     refs=[("A. Almenningen, O. Bastiansen, V. Ewing, K. Hedberg, M. Trætteberg, The molecular structure of " *
            "disiloxane, (SiH3)2O, Acta Chem. Scand. 17, 2455 (1963)", "10.3891/acta.chem.scand.17-2455")],
     note="Electron diffraction. The wide Si–O–Si angle (144°) is a thermal average of a floppy bend; the SiH3 " *
          "orientation (one H of each group towards the C2 axis) is the authors' tentative choice."),
    (name="Trisilane", formula="Si3H8", category="Inorganic: hydrides", use="deposition precursor", pg="C2v",
     geometry=trisilane("Si", 2.332, 1.483, 1.486, 110.2, 112.9, 107.8),
     refs=[("A. Haaland, K. Rypdal, H. Stüger, H. V. Volden, Molecular structures and conformational composition of " *
            "trisilane and tetrasilane by gas-phase electron diffraction, Acta Chem. Scand. 48, 46 (1994)",
            "10.3891/acta.chem.scand.48-0046")],
     note="Electron-diffraction r_a (Table 1); the central H–Si–H angle (107.8°) was fixed at an ab initio value " *
          "by the authors."),
    # Unsaturated organics, cyanides, interhalogens
    (name="Diacetylene", formula="C4H2", category="Organic: hydrocarbons", pg="D∞h",
     geometry=linear(("H", "C", "C", "C", "C", "H"), (1.0935, 1.2176, 1.3837, 1.2176, 1.0935)),
     refs=[("M. Tanimoto, K. Kuchitsu, Y. Morino, Molecular structure of diacetylene as studied by gas electron " *
            "diffraction, Bull. Chem. Soc. Jpn. 44, 386 (1971)", "10.1246/bcsj.44.386")],
     note="Electron-diffraction r_g (thermal average)."),
    (name="Cyclopropene", formula="C3H4", category="Organic: hydrocarbons", pg="C2v",
     geometry=cyclopropene(1.509, 1.2959, 1.088, 1.072, 149.85, 114.57),
     refs=[("W. M. Stigliani, V. W. Laurie, J. C. Li, Microwave spectrum and structure of cyclopropene, J. Chem. " *
            "Phys. 62, 1890 (1975)", "10.1063/1.430674")],
     note="Microwave structure (13C species). The ring is placed from the two C–C distances, giving a C–C–C " *
          "angle of 50.86° (authors: 50.84 ± 0.05°)."),
    (name="Fluoroacetylene", formula="HCCF", category="Organic: halocarbons", pg="C∞v",
     geometry=linear(("H", "C", "C", "F"), (1.053, 1.198, 1.279)), refs=[TYLER_SHERIDAN],
     note=TYLER_SHERIDAN_NOTE * " F is placed from the centre-of-mass condition."),
    (name="Chloroacetylene", formula="HCCCl", category="Organic: halocarbons", pg="C∞v",
     geometry=linear(("H", "C", "C", "Cl"), (1.0550, 1.2033, 1.6368)),
     refs=[("H. Jones, M. Takami, J. Sheridan, Microwave spectrum of chloroacetylene in ground and excited " *
            "vibrational states, Z. Naturforsch. A 33, 156 (1978)", "10.1515/zna-1978-0209")],
     note="Substitution structure r_s, average over four isotopologues (consistent to 0.0005 Å)."),
    (name="Vinyl chloride", formula="C2H3Cl", category="Organic: halocarbons", pg="Cs",
     geometry=planar_ab("C", "C", 1.332, [("Cl", 1, 1.726, 122.0, 1), ("H", 1, 1.08, 124.0, -1),
                                         ("H", 2, 1.08, 121.0, 1), ("H", 2, 1.08, 120.0, -1)]),
     refs=[("D. Kivelson, E. B. Wilson, D. R. Lide, Microwave spectrum, structure, dipole moment, and nuclear " *
            "quadrupole effects in vinyl chloride, J. Chem. Phys. 32, 205 (1960)", "10.1063/1.1700901")],
     note="Microwave structure from nine isotopologues, planarity assumed; the authors give one average C–H " *
          "(1.08 Å), used for all three, and angles to 1°."),
    (name="1,1-Difluoroethylene", formula="H2CCF2", category="Organic: halocarbons", pg="C2v",
     geometry=planar_ab("C", "C", 1.315, [("F", 1, 1.323, 180 - 109.1 / 2, 1), ("F", 1, 1.323, 180 - 109.1 / 2, -1),
                                         ("H", 2, 1.079, 180 - 121.8 / 2, 1), ("H", 2, 1.079, 180 - 121.8 / 2, -1)]),
     refs=[("V. W. Laurie, D. T. Pence, Microwave spectra and structures of difluoroethylenes, J. Chem. Phys. 38, " *
            "2693 (1963)", "10.1063/1.1733575")],
     note="Effective microwave structure (vinylidene fluoride), no uncertainties given."),
    (name="cis-1,2-Difluoroethylene", formula="CHFCHF", category="Organic: halocarbons", pg="C2v",
     geometry=planar_ab("C", "C", 1.324, [("F", 1, 1.335, 122.1, 1), ("F", 2, 1.335, 122.1, 1),
                                         ("H", 1, 1.089, 124.0, -1), ("H", 2, 1.089, 124.0, -1)]),
     refs=[("V. W. Laurie, D. T. Pence, Microwave spectra and structures of difluoroethylenes, J. Chem. Phys. 38, " *
            "2693 (1963)", "10.1063/1.1733575")],
     note="Effective microwave structure, no uncertainties given."),
    (name="Methanimine", formula="CH2NH", category="Organic: nitrogen and sulfur compounds", pg="Cs",
     geometry=planar_ab("C", "N", 1.273, [("H", 2, 1.021, 110.4, 1), ("H", 1, 1.09, 125.1, 1),
                                         ("H", 1, 1.09, 360 - 125.1 - 117.0, -1)]),
     refs=[("R. Pearson, F. J. Lovas, Microwave spectrum and molecular structure of methylenimine (CH2NH), J. Chem. " *
            "Phys. 66, 4149 (1977)", "10.1063/1.434490")],
     note="Pseudo-Kraitchman structure with the two C–H distances constrained equal by the authors; the C–H " *
          "cis to the N–H at N=C–H 125.1°, H–C–H 117.0°."),
    (name="Cyanogen", formula="C2N2", category="Organic: nitrogen and sulfur compounds", pg="D∞h",
     geometry=linear(("N", "C", "C", "N"), (1.154, 1.389, 1.154)),
     refs=[("A. G. Maki, High-resolution infrared spectra of cyanogen and cyanogen-15N2, J. Chem. Phys. 43, 3193 " *
            "(1965)", "10.1063/1.1697293")],
     note="Ground-state r_0 of low precision (3σ: C–C ± 0.030 Å, C≡N ± 0.017 Å)."),
    (name="Cyanoacetylene", formula="HC3N", category="Organic: nitrogen and sulfur compounds", pg="C∞v",
     geometry=linear(("H", "C", "C", "C", "N"), (1.058, 1.205, 1.378, 1.159)), refs=[TYLER_SHERIDAN],
     note=TYLER_SHERIDAN_NOTE),
    (name="Cyanogen fluoride", formula="FCN", category="Organic: nitrogen and sulfur compounds", pg="C∞v",
     geometry=linear(("F", "C", "N"), (1.262, 1.159)), refs=[TYLER_SHERIDAN],
     note=TYLER_SHERIDAN_NOTE * " F is placed from the centre-of-mass condition."),
    (name="Cyanogen chloride", formula="ClCN", category="Organic: nitrogen and sulfur compounds", pg="C∞v",
     geometry=linear(("Cl", "C", "N"), (1.631, 1.159)), refs=[TYLER_SHERIDAN], note=TYLER_SHERIDAN_NOTE),
    (name="Cyanogen bromide", formula="BrCN", category="Organic: nitrogen and sulfur compounds", pg="C∞v",
     geometry=linear(("Br", "C", "N"), (1.789, 1.158)), refs=[TYLER_SHERIDAN], note=TYLER_SHERIDAN_NOTE),
    (name="Cyanogen iodide", formula="ICN", category="Organic: nitrogen and sulfur compounds", pg="C∞v",
     geometry=linear(("I", "C", "N"), (1.994, 1.159)), refs=[TYLER_SHERIDAN],
     note=TYLER_SHERIDAN_NOTE * " I is placed from the centre-of-mass condition."),
    (name="Carbon suboxide", formula="C3O2", category="Inorganic: oxides and sulfides", pg="D∞h",
     geometry=linear(("O", "C", "C", "C", "O"), (1.1632, 1.2894, 1.2894, 1.1632)),
     refs=[("M. Tanimoto, K. Kuchitsu, Y. Morino, Molecular structure and the effect of large-amplitude vibration " *
            "of carbon suboxide as studied by gas electron diffraction, Bull. Chem. Soc. Jpn. 43, 2776 (1970)",
            "10.1246/bcsj.43.2776")],
     note="Electron-diffraction r_g. Quasi-linear: the C=C=C angle bends with a very large amplitude (about " *
          "68 cm⁻¹), so the linear chain is a model of the average."),
    (name="Chlorine dioxide", formula="ClO2", category="Inorganic: oxides and sulfides", pg="C2v",
     geometry=zmatrix([("Cl",), ("O", 1, 1.469839), ("O", 1, 1.469839, 2, 117.4033)]),
     refs=[("K. Miyazaki, M. Tanoura, K. Tanaka, T. Tanaka, Microwave spectrum of chlorine dioxide in excited " *
            "vibrational states, J. Mol. Spectrosc. 116, 435 (1986)", "10.1016/0022-2852(86)90138-4")],
     note="Equilibrium r_e, read in the authors' abstract (CiNii). " * odd_electron_note()),
    (name="Bromine monochloride", formula="BrCl", category="Inorganic: nonmetal halides", pg="C∞v",
     geometry=linear(("Br", "Cl"), (2.13605328,)),
     refs=[("H. Uehara, T. Konno, Y. Ozaki, K. Horiai, K. Nakagawa, J. W. C. Johns, The high-resolution infrared " *
            "spectrum of BrCl, Can. J. Phys. 72, 1145 (1994)", "10.1139/p94-148")],
     note="Equilibrium r_e."),
    (name="Bromine trifluoride", formula="BrF3", category="Inorganic: nonmetal halides", pg="C2v",
     geometry=zmatrix([("Br",), ("F", 1, 1.721), ("F", 1, 1.810, 2, 86.21), ("F", 1, 1.810, 2, 86.21, 3, 180.0)]),
     refs=[("D. W. Magnuson, Microwave spectrum and molecular structure of bromine trifluoride, J. Chem. Phys. 27, " *
            "223 (1957)", "10.1063/1.1743675")],
     note="Planar T shape like ClF3; effective ground-state structure from the moments of inertia. The " *
          "F–Br–F′ angle (86° 12.6′) is taken between the short bond and each long bond."),
    (name="Iodine heptafluoride", formula="IF7", category="Inorganic: nonmetal halides", pg="D5h",
     geometry=pentagonal_bipyramid("I", "F", 1.786, 1.858),
     refs=[("W. J. Adams, H. B. Thompson, L. S. Bartell, Structure, pseudorotation, and vibrational mode coupling " *
            "in IF7: an electron diffraction study, J. Chem. Phys. 53, 4040 (1970)", "10.1063/1.1673876")],
     note="Electron diffraction, idealised D5h pentagonal bipyramid. The molecule pseudorotates: on average the " *
          "ring puckers by 7.5° and the axial bonds bend by 4.5°."),
    # Halomethanes and fluorocarbon etch gases
    (name="Bromomethane", formula="CH3Br", category="Organic: halocarbons", pg="C3v",
     geometry=xy3z("C", "H", "Br", 1.086, dm(111, 10), 1.933),
     refs=[("J. L. Duncan, The equilibrium geometries of the methyl halides, J. Mol. Struct. 6, 447 (1970)",
            "10.1016/0022-2860(70)90027-X")],
     note="Equilibrium r_e, read in the authors' abstract."),
    (name="Trifluoroiodomethane", formula="CF3I", category="Organic: halocarbons", use="etch gas", pg="C3v",
     geometry=xy3z("C", "F", "I", 1.3285, 108.42, 2.1438),
     refs=[("A. P. Cox, G. Duxbury, J. A. Hardy, Y. Kawashima, Microwave spectra of CF3Br and CF3I. Structures and " *
            "dipole moments, J. Chem. Soc., Faraday Trans. 2 76, 339 (1980)", "10.1039/F29807600339")],
     note="Zero-point average structure r_z."),
    (name="Bromotrifluoromethane", formula="CF3Br", category="Organic: halocarbons", pg="C3v",
     geometry=xy3z("C", "F", "Br", 1.3265, 108.81, 1.9234),
     refs=[("A. P. Cox, G. Duxbury, J. A. Hardy, Y. Kawashima, Microwave spectra of CF3Br and CF3I. Structures and " *
            "dipole moments, J. Chem. Soc., Faraday Trans. 2 76, 339 (1980)", "10.1039/F29807600339")],
     note="Zero-point average structure r_z."),
    (name="Chlorotrifluoromethane", formula="CF3Cl", category="Organic: halocarbons", pg="C3v",
     geometry=xy3z("C", "F", "Cl", 1.328, 108.6, 1.751),
     refs=[("L. S. Bartell, L. O. Brockway, Structures of CF3Cl and CH3Cl, J. Chem. Phys. 23, 1860 (1955)",
            "10.1063/1.1740593")],
     note="Electron diffraction, most probable distances; no uncertainties given in the abstract."),
    (name="Dichlorodifluoromethane", formula="CF2Cl2", category="Organic: halocarbons", pg="C2v",
     geometry=xy2z2("C", "Cl", "F", 1.744, dm(112, 33), 1.345, dm(106, 14)),
     refs=[("H. Takeo, C. Matsumura, The microwave spectrum of dichlorodifluoromethane, Bull. Chem. Soc. Jpn. 50, " *
            "636 (1977)", "10.1246/bcsj.50.636")],
     note="Microwave structure from three isotopologues; no uncertainties given."),
    (name="Chlorodifluoromethane", formula="CHClF2", category="Organic: halocarbons", pg="Cs",
     geometry=xzy2w("C", "Cl", "F", "H", 1.74, 1.35, 1.09, dm(110, 30), 107.0, 107.0),
     refs=[("E. L. Beeson, T. L. Weatherly, Q. Williams, Molecular constants of chlorodifluoromethane from " *
            "microwave spectrum analysis, J. Chem. Phys. 37, 2926 (1962)", "10.1063/1.1733119")],
     note="Low-precision microwave structure (two isotopologues); the C–H distance (1.09 Å) was assumed by the " *
          "authors."),
    (name="Fluoroform", formula="CHF3", category="Organic: halocarbons", use="etch gas", pg="C3v",
     geometry=xy3z("C", "F", "H", 1.332, dm(108, 48), 1.098),
     refs=[("S. N. Ghosh, R. Trambarulo, W. Gordy, Microwave spectra and molecular structures of fluoroform, " *
            "chloroform, and methyl chloroform, J. Chem. Phys. 20, 605 (1952)", "10.1063/1.1700501")],
     note="Microwave structure from several isotopic species; no uncertainties given."),
    (name="Difluoromethane", formula="CH2F2", category="Organic: halocarbons", use="etch gas", pg="C2v",
     geometry=xy2z2("C", "F", "H", 1.3508, 108.49, 1.084, 112.8),
     refs=[("E. Hirota, Anharmonic potential function and equilibrium structure of methylene fluoride, J. Mol. " *
            "Spectrosc. 71, 145 (1978)", "10.1016/0022-2852(78)90079-6")],
     note="Equilibrium r_e from five isotopologues, read in the authors' abstract."),
    (name="Dibromomethane", formula="CH2Br2", category="Organic: halocarbons", pg="C2v",
     geometry=xy2z2("C", "Br", "H", 1.930, 112.5, 1.079, 113.6),
     refs=[("D. Chadwick, D. J. Millen, Microwave spectrum of dibromomethane. Part 1. Molecular structure, Trans. " *
            "Faraday Soc. 67, 1539 (1971)", "10.1039/TF9716701539")],
     note="Ground-state r_0, the authors' structure (I) with all parameters fitted (C–H ± 0.01 Å, H–C–H ± 1.5°)."),
    (name="Bromoform", formula="CHBr3", category="Organic: halocarbons", pg="C3v",
     geometry=xy3z("C", "Br", "H", 1.924, 111.7, 1.111),
     refs=[("K. Tamagawa, M. Kimura, Molecular structure of bromoform as determined by a joint analysis of electron " *
            "diffraction and microwave data, Bull. Chem. Soc. Jpn. 52, 2747 (1979)", "10.1246/bcsj.52.2747")],
     note="Electron diffraction with microwave data: r_g distances (C–H ± 0.046 Å), average angles."),
    (name="Difluorocarbene", formula="CF2", category="Organic: radicals", pg="C2v",
     geometry=zmatrix([("C",), ("F", 1, 1.300), ("F", 1, 1.300, 2, 104.94)]),
     refs=[("C. W. Mathews, The absorption spectrum of CF2, Can. J. Phys. 45, 2355 (1967)", "10.1139/p67-188")],
     note="Ground-state r_0 from the electronic spectrum; no uncertainties given. " * SINGLET_NOTE),
    (name="Dichlorocarbene", formula="CCl2", category="Organic: radicals", pg="C2v",
     geometry=zmatrix([("C",), ("Cl", 1, 1.7157), ("Cl", 1, 1.7157, 2, 109.2)]),
     refs=[("M. Fujitake, E. Hirota, The millimeter- and submillimeter-wave spectrum of dichlorocarbene CCl2, " *
            "J. Chem. Phys. 91, 3426 (1989)", "10.1063/1.456916")],
     note="Ground-state r_0. " * SINGLET_NOTE),
    # Heavier hydrides, small inorganic molecules, Si/Ge halides, high-temperature diatomics
    (name="Hydrogen selenide", formula="H2Se", category="Inorganic: hydrides", pg="C2v",
     geometry=zmatrix([("Se",), ("H", 1, 1.460), ("H", 1, 1.460, 2, dm(90, 34))]),
     refs=[("R. A. Hill, T. H. Edwards, Vibrational analysis and isotope effects in hydrogen selenide, J. Chem. " *
            "Phys. 42, 1391 (1965)", "10.1063/1.1696127")],
     note="Equilibrium structure calculated by the authors from experimental vibration–rotation constants; no " *
          "uncertainties given."),
    (name="Disulfane", formula="H2S2", category="Inorganic: hydrides", pg="C2",
     geometry=zmatrix([("S",), ("S", 1, 2.0564), ("H", 1, 1.3421, 2, 97.88), ("H", 2, 1.3421, 1, 97.88, 3, 90.34)]),
     refs=[("J. Hahn, P. Schmidt, K. Reinartz, J. Behrend, G. Winnewisser, K. M. T. Yamada, Synthesis and molecular " *
            "structure of disulfane, Z. Naturforsch. B 46, 1338 (1991)", "10.1515/znb-1991-1011")],
     note="The authors' \"partially corrected equilibrium structure\"; dihedral H–S–S–H 90.3(2)°."),
    (name="Methylidynephosphine", formula="HCP", category="Inorganic: hydrides", pg="C∞v",
     geometry=linear(("H", "C", "P"), (1.0666, 1.54020)),
     refs=[("J. Lavigne, A. Cabana, The ν3 vibration–rotation band of DCP and the equilibrium structure of " *
            "methinophosphide, Can. J. Phys. 60, 304 (1982)", "10.1139/p82-039")],
     note="Equilibrium r_e. The phosphorus analogue of HCN (phosphaethyne)."),
    (name="Hydrogen isocyanide", formula="HNC", category="Organic: nitrogen and sulfur compounds", pg="C∞v",
     geometry=linear(("H", "N", "C"), (0.987, 1.171)),
     refs=[("R. A. Creswell, E. F. Pearson, M. Winnewisser, G. Winnewisser, Detection of the millimeter wave " *
            "spectrum of hydrogen isocyanide, HNC, Z. Naturforsch. A 31, 221 (1976)", "10.1515/zna-1976-3-401")],
     note="Substitution-type structure (H and C by substitution, N from the centre of mass); ± 0.003–0.004 Å."),
    (name="Nitroxyl", formula="HNO", category="Inorganic: oxoacids", pg="Cs",
     geometry=zmatrix([("N",), ("O", 1, 1.2058), ("H", 1, 1.0628, 2, 109.09)]),
     refs=[("E. Hirota, A simultaneous determination of third-order vibrational anharmonicity constants and " *
            "equilibrium structure. Molecular structure of non-linear HXY-type molecules, Nippon Kagaku Kaishi " *
            "1986, 1438", "10.1246/nikkashi.1986.1438")],
     note="Equilibrium r_e (3σ: N–H ± 0.0025 Å), read in the author's abstract (CiNii)."),
    (name="Hypochlorous acid", formula="HOCl", category="Inorganic: oxoacids", pg="Cs",
     geometry=zmatrix([("O",), ("Cl", 1, 1.6891), ("H", 1, 0.9654, 2, 103.21)]),
     refs=[("E. Hirota, A simultaneous determination of third-order vibrational anharmonicity constants and " *
            "equilibrium structure. Molecular structure of non-linear HXY-type molecules, Nippon Kagaku Kaishi " *
            "1986, 1438", "10.1246/nikkashi.1986.1438")],
     note="Equilibrium r_e (3σ: O–H ± 0.0035 Å), read in the author's abstract (CiNii)."),
    (name="Difluorosilane", formula="SiH2F2", category="Inorganic: nonmetal halides", pg="C2v",
     geometry=xy2z2("Si", "F", "H", 1.5767, dm(107, 56), 1.471, dm(112, 1)),
     refs=[("V. W. Laurie, Microwave spectrum, dipole moment, and structure of difluorosilane, J. Chem. Phys. 26, " *
            "1359 (1957)", "10.1063/1.1743546")],
     note="Ground-state microwave structure (Si–H ± 0.007 Å)."),
    (name="Trifluorosilane", formula="SiHF3", category="Inorganic: nonmetal halides", pg="C3v",
     geometry=xy3z("Si", "F", "H", 1.561, dm(108, 6), 1.55),
     refs=[("J. Sheridan, W. Gordy, The microwave spectra and molecular structures of trifluorosilane derivatives, " *
            "J. Chem. Phys. 19, 965 (1951)", "10.1063/1.1748418")],
     note="Microwave structure in which the authors assumed Si–H = 1.55 ± 0.05 Å; Si–F and F–Si–F follow from it."),
    (name="Silicon difluoride", formula="SiF2", category="Inorganic: nonmetal halides", pg="C2v",
     geometry=zmatrix([("Si",), ("F", 1, 1.5947), ("F", 1, 1.5947, 2, dm(100, 53))]),
     refs=[("V. M. Khanna, R. Hauge, R. F. Curl, J. L. Margrave, Infrared spectrum, force constants, and " *
            "thermodynamic functions of SiF2, J. Chem. Phys. 47, 5031 (1967)", "10.1063/1.1701755")],
     note="Ground-state average structure ⟨r⟩ (r_z type); no uncertainties given. A silylene: " * SINGLET_NOTE),
    (name="Germanium difluoride", formula="GeF2", category="Inorganic: nonmetal halides", pg="C2v",
     geometry=zmatrix([("Ge",), ("F", 1, 1.736), ("F", 1, 1.736, 2, 97.4)]),
     refs=[("H. Takeo, R. F. Curl, P. W. Wilson, Microwave spectrum of GeF2, Ohio State University Symposium on " *
            "Molecular Spectroscopy, 1970, http://hdl.handle.net/1811/8588", "")],
     note="Ground-state structure from the authors' symposium abstract; no uncertainties given."),
    (name="Germanium dichloride", formula="GeCl2", category="Inorganic: nonmetal halides", pg="C2v",
     geometry=zmatrix([("Ge",), ("Cl", 1, 2.169452), ("Cl", 1, 2.169452, 2, 99.8825)]),
     refs=[("M. J. Tsuchiya, H. Honjou, K. Tanaka, T. Tanaka, Millimeter-wave spectrum of germanium dichloride GeCl2. " *
            "Equilibrium structure and anharmonic force field, J. Mol. Struct. 352–353, 407 (1995)",
            "10.1016/0022-2860(95)08830-O")],
     note="Equilibrium r_e, read in the authors' abstract (CiNii)."),
    (name="Boron monofluoride", formula="BF", category="Inorganic: nonmetal halides", pg="C∞v",
     geometry=linear(("B", "F"), (1.2625,)),
     refs=[("F. J. Lovas, D. R. Johnson, Microwave spectrum of BF, J. Chem. Phys. 55, 41 (1971)", "10.1063/1.1675537")],
     note="Equilibrium r_e (11BF); a high-temperature species."),
    (name="Aluminium monofluoride", formula="AlF", category="Inorganic: metal halides", pg="C∞v",
     geometry=linear(("Al", "F"), (1.6543601,)),
     refs=[("J. Hoeft, F. J. Lovas, E. Tiemann, T. Törring, Microwave absorption spectra of AlF, GaF, InF, and TlF, " *
            "Z. Naturforsch. A 25, 1029 (1970)", "10.1515/zna-1970-0706")],
     note="Equilibrium r_e (Table 13); a high-temperature species."),
    (name="Aluminium monochloride", formula="AlCl", category="Inorganic: metal halides", pg="C∞v",
     geometry=linear(("Al", "Cl"), (2.13011,)),
     refs=[("F. C. Wyse, W. Gordy, Millimeter wave rotational spectra of AlCl, AlBr, and AlI, J. Chem. Phys. 56, " *
            "2130 (1972)", "10.1063/1.1677509")],
     note="Equilibrium r_e; a high-temperature species."),
    (name="Gallium monochloride", formula="GaCl", category="Inorganic: metal halides", pg="C∞v",
     geometry=linear(("Ga", "Cl"), (2.201681,)),
     refs=[("E. Tiemann, M. Grasshoff, J. Hoeft, Hyperfeinstruktur im Rotationsspektrum des GaCl, Z. Naturforsch. A " *
            "27, 753 (1972)", "10.1515/zna-1972-0505")],
     note="Equilibrium r_e of 69Ga35Cl (Tab. 5); the gallium carrier in GaCl-based (HVPE) growth."),
    (name="Silicon monosulfide", formula="SiS", category="Inorganic: oxides and sulfides", pg="C∞v",
     geometry=linear(("Si", "S"), (1.9293141,)),
     refs=[("J. Hoeft, F. J. Lovas, T. Törring, Elektrisches Dipolmoment und Mikrowellenrotationsspektrum von SiS, " *
            "Z. Naturforsch. A 24, 1422 (1969)", "10.1515/zna-1969-0932")],
     note="Equilibrium r_e (Table 4)."),
    (name="Germanium monoxide", formula="GeO", category="Inorganic: oxides and sulfides", pg="C∞v",
     geometry=linear(("Ge", "O"), (1.624647,)),
     refs=[("T. Törring, Das Mikrowellenrotationsspektrum des GeO, Z. Naturforsch. A 21, 287 (1966)",
            "10.1515/zna-1966-0319")],
     note="Equilibrium r_e (Tab. 2), 1Σ ground state."),
    (name="Phosphorus mononitride", formula="PN", category="Inorganic: nitrides", pg="C∞v",
     geometry=linear(("P", "N"), (1.49085,)),
     refs=[("J. Hoeft, E. Tiemann, T. Törring, Rotationsspektrum des PN, Z. Naturforsch. A 27, 703 (1972)",
            "10.1515/zna-1972-0424")],
     note="Equilibrium r_e (Tab. 2)."),
    (name="Silicon monofluoride", formula="SiF", category="Inorganic: radicals", pg="C∞v",
     geometry=linear(("Si", "F"), (1.601018,)),
     refs=[("K. Tanaka, Y. Akiyama, T. Tanaka, Infrared diode laser spectroscopy of the SiF radical, J. Mol. " *
            "Spectrosc. 137, 55 (1989)", "10.1016/0022-2852(89)90268-3")],
     note="Equilibrium r_e of the X 2Π ground state, read in the authors' abstract (CiNii). " * spin_note("doublet")),
    (name="Phosphino radical", formula="PH2", category="Inorganic: radicals", pg="C2v",
     geometry=zmatrix([("P",), ("H", 1, 1.43365), ("H", 1, 1.43365, 2, 91.622)]),
     refs=[("T. Hirao, S. Hayakashi, S. Yamamoto, S. Saito, Microwave spectrum of the PD2 radical in the 2B1 ground " *
            "electronic state, J. Mol. Spectrosc. 187, 153 (1998)", "10.1006/jmsp.1997.7457")],
     note="Zero-point average structure r_z of the X̃ 2B1 ground state, read in the authors' abstract (CiNii). " *
          spin_note("doublet")),
    # Carbon nanostructures
    (name="Buckminsterfullerene", formula="C60", category="Carbon: fullerenes", pg="Ih",
     geometry=fullerene_c60(1.458, 1.401),
     refs=[("K. Hedberg, L. Hedberg, D. S. Bethune, C. A. Brown, H. C. Dorn, R. D. Johnson, M. de Vries, Bond " *
            "lengths in free molecules of buckminsterfullerene, C60, from gas-phase electron diffraction, Science 254, " *
            "410 (1991)", "10.1126/science.254.5030.410")],
     note="Electron-diffraction r_g (thermal average at about 1000 K): 1.458 Å in the pentagons (5:6 bonds), " *
          "1.401 Å between hexagons (6:6 bonds)."),
]

"""
Coordinates centred and turned to their principal axes (largest extent along
x, smallest along z, so planar molecules lie in the xy plane).
"""
function orient(pos)
    P = reduce(hcat, pos)
    P = P .- sum(P; dims=2) / size(P, 2)
    F = eigen(Symmetric(P * P'))
    R = F.vectors[:, sortperm(F.values; rev=true)]
    for k in 1:3                                           # fixed signs, right-handed
        v = R[:, k]
        i = findfirst(x -> abs(x) > 1e-8, v)
        i !== nothing && v[i] < 0 && (R[:, k] = -v)
    end
    det(R) < 0 && (R[:, 3] = -R[:, 3])
    Q = R' * P
    return [Q[:, k] .+ 0.0 for k in 1:size(Q, 2)]
end
