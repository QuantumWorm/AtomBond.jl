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

spin_note(state) = "The ground state is a spin $state (open shell): use a spin-polarized calculation."

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
    (name="Carbon monoxide", formula="CO", category="Inorganic: oxides", pg="C∞v",
     geometry=linear(("C", "O"), (1.1283,)), refs=[HUBER_HERZBERG]),
    (name="Nitric oxide", formula="NO", category="Inorganic: oxides", pg="C∞v",
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
    (name="Carbon dioxide", formula="CO2", category="Inorganic: oxides", pg="D∞h",
     geometry=linear(("O", "C", "O"), (1.1621, 1.1621)), refs=[HERZBERG_1966]),
    (name="Hydrogen sulfide", formula="H2S", category="Inorganic: hydrides", pg="C2v",
     geometry=zmatrix([("S",), ("H", 1, 1.3356), ("H", 1, 1.3356, 2, 92.11)]),
     refs=[("R. L. Cook, F. C. De Lucia, P. Helminger, Molecular force field and structure of hydrogen sulfide: " *
            "recent microwave results, J. Mol. Struct. 28, 237 (1975)", "10.1016/0022-2860(75)80094-9")]),
    (name="Nitrous oxide", formula="N2O", category="Inorganic: oxides", use="ALD/CVD co-reactant (oxidant, nitrogen source)", pg="C∞v",
     geometry=linear(("N", "N", "O"), (1.1282, 1.1842)), refs=[HERZBERG_1966]),
    (name="Sulfur dioxide", formula="SO2", category="Inorganic: oxides", pg="C2v",
     geometry=zmatrix([("S",), ("O", 1, 1.4321), ("O", 1, 1.4321, 2, 119.5)]), refs=[HERZBERG_1966]),
    (name="Ozone", formula="O3", category="Inorganic: elements", use="ALD/CVD co-reactant (oxidant, nitrogen source)", pg="C2v",
     geometry=zmatrix([("O",), ("O", 1, 1.278), ("O", 1, 1.278, 2, 116.8)]), refs=[HERZBERG_1966]),
    (name="Hydrogen peroxide", formula="H2O2", category="Inorganic: hydrides", pg="C2",
     geometry=zmatrix([("O",), ("O", 1, 1.475), ("H", 1, 0.950, 2, 94.8), ("H", 2, 0.950, 1, 94.8, 3, 119.8)]),
     refs=[("R. L. Redington, W. B. Olson, P. C. Cross, Studies of hydrogen peroxide: the infrared spectrum and " *
            "the internal rotation problem, J. Chem. Phys. 36, 1311 (1962)", "10.1063/1.1732733")]),
    (name="Hydrogen cyanide", formula="HCN", category="Organic: nitrogen and sulfur compounds", pg="C∞v",
     geometry=linear(("H", "C", "N"), (1.064, 1.156)), refs=[HERZBERG_1966]),
    (name="Carbonyl sulfide", formula="OCS", category="Inorganic: oxides", pg="C∞v",
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
