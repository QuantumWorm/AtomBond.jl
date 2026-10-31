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

# --- references --------------------------------------------------------------------------

const HUBER_HERZBERG = ("K. P. Huber, G. Herzberg, Molecular Spectra and Molecular Structure IV. Constants of " *
                        "Diatomic Molecules, Van Nostrand Reinhold (1979)", "10.1007/978-1-4757-0961-2")
const HERZBERG_1966 = ("G. Herzberg, Molecular Spectra and Molecular Structure III. Electronic Spectra and " *
                       "Electronic Structure of Polyatomic Molecules, Van Nostrand (1966)", "")
const KUCHITSU = ("K. Kuchitsu (ed.), Structure of Free Polyatomic Molecules - Basic Data, Springer (1998)",
                  "10.1007/978-3-642-45748-7")
const GURVICH = ("L. V. Gurvich, I. V. Veyts, C. B. Alcock, Thermodynamic Properties of Individual Substances, " *
                 "4th ed., Hemisphere (1989)", "")

spin_note(state) = "The ground state is a spin $state (open shell): use a spin-polarized calculation."

# name, formula, category, point group, geometry (symbols, positions), references
# (text, DOI), optional note. The values are those of the references.
const MOLECULES = [
    # Diatomic molecules (equilibrium bond lengths r_e)
    (name="Hydrogen", formula="H2", category="Diatomic molecules", pg="D∞h",
     geometry=linear(("H", "H"), (0.7414,)), refs=[HUBER_HERZBERG]),
    (name="Nitrogen", formula="N2", category="Diatomic molecules", pg="D∞h",
     geometry=linear(("N", "N"), (1.0977,)), refs=[HUBER_HERZBERG]),
    (name="Oxygen", formula="O2", category="Diatomic molecules", pg="D∞h",
     geometry=linear(("O", "O"), (1.2075,)), refs=[HUBER_HERZBERG], note=spin_note("triplet")),
    (name="Carbon monoxide", formula="CO", category="Diatomic molecules", pg="C∞v",
     geometry=linear(("C", "O"), (1.1283,)), refs=[HUBER_HERZBERG]),
    (name="Nitric oxide", formula="NO", category="Diatomic molecules", pg="C∞v",
     geometry=linear(("N", "O"), (1.1508,)), refs=[HUBER_HERZBERG], note=spin_note("doublet")),
    (name="Hydrogen fluoride", formula="HF", category="Diatomic molecules", pg="C∞v",
     geometry=linear(("H", "F"), (0.9168,)),
     refs=[("R. J. Le Roy, Improved parameterization for combined isotopomer analysis of diatomic spectra and its " *
            "application to HF and DF, J. Mol. Spectrosc. 194, 189 (1999)", "10.1006/jmsp.1998.7786")]),
    (name="Hydrogen chloride", formula="HCl", category="Diatomic molecules", pg="C∞v",
     geometry=linear(("H", "Cl"), (1.2746,)), refs=[HUBER_HERZBERG]),
    # Small molecules
    (name="Water", formula="H2O", category="Small molecules", pg="C2v",
     geometry=zmatrix([("O",), ("H", 1, 0.9578), ("H", 1, 0.9578, 2, 104.48)]),
     refs=[("A. R. Hoy, P. R. Bunker, A precise solution of the rotation bending Schrödinger equation for a " *
            "triatomic molecule with application to the water molecule, J. Mol. Spectrosc. 74, 1 (1979)",
            "10.1016/0022-2852(79)90019-5"), HERZBERG_1966]),
    (name="Ammonia", formula="NH3", category="Small molecules", pg="C3v",
     geometry=pyramid("N", "H", 1.0124, 106.67), refs=[HERZBERG_1966]),
    (name="Carbon dioxide", formula="CO2", category="Small molecules", pg="D∞h",
     geometry=linear(("O", "C", "O"), (1.1621, 1.1621)), refs=[HERZBERG_1966]),
    (name="Hydrogen sulfide", formula="H2S", category="Small molecules", pg="C2v",
     geometry=zmatrix([("S",), ("H", 1, 1.3356), ("H", 1, 1.3356, 2, 92.11)]),
     refs=[("R. L. Cook, F. C. De Lucia, P. Helminger, Molecular force field and structure of hydrogen sulfide: " *
            "recent microwave results, J. Mol. Struct. 28, 237 (1975)", "10.1016/0022-2860(75)80094-9")]),
    (name="Nitrous oxide", formula="N2O", category="Small molecules", pg="C∞v",
     geometry=linear(("N", "N", "O"), (1.1282, 1.1842)), refs=[HERZBERG_1966]),
    (name="Sulfur dioxide", formula="SO2", category="Small molecules", pg="C2v",
     geometry=zmatrix([("S",), ("O", 1, 1.4321), ("O", 1, 1.4321, 2, 119.5)]), refs=[HERZBERG_1966]),
    (name="Ozone", formula="O3", category="Small molecules", pg="C2v",
     geometry=zmatrix([("O",), ("O", 1, 1.278), ("O", 1, 1.278, 2, 116.8)]), refs=[HERZBERG_1966]),
    (name="Hydrogen peroxide", formula="H2O2", category="Small molecules", pg="C2",
     geometry=zmatrix([("O",), ("O", 1, 1.475), ("H", 1, 0.950, 2, 94.8), ("H", 2, 0.950, 1, 94.8, 3, 119.8)]),
     refs=[("R. L. Redington, W. B. Olson, P. C. Cross, Studies of hydrogen peroxide: the infrared spectrum and " *
            "the internal rotation problem, J. Chem. Phys. 36, 1311 (1962)", "10.1063/1.1732733")]),
    (name="Hydrogen cyanide", formula="HCN", category="Small molecules", pg="C∞v",
     geometry=linear(("H", "C", "N"), (1.064, 1.156)), refs=[HERZBERG_1966]),
    # Organic molecules
    (name="Methane", formula="CH4", category="Organic molecules", pg="Td",
     geometry=tetrahedron("C", "H", 1.087),
     refs=[("E. Hirota, Anharmonic potential function and equilibrium structure of methane, J. Mol. Spectrosc. " *
            "77, 213 (1979)", "10.1016/0022-2852(79)90103-6"), HERZBERG_1966]),
    (name="Ethane", formula="C2H6", category="Organic molecules", pg="D3d",
     geometry=staggered("C", "H", 1.536, 1.091, 110.91), refs=[HERZBERG_1966, KUCHITSU]),
    (name="Ethylene", formula="C2H4", category="Organic molecules", pg="D2h",
     geometry=planar_x2y4("C", "H", 1.339, 1.086, 117.6), refs=[HERZBERG_1966]),
    (name="Acetylene", formula="C2H2", category="Organic molecules", pg="D∞h",
     geometry=linear(("H", "C", "C", "H"), (1.0631, 1.2026, 1.0631)), refs=[HERZBERG_1966, KUCHITSU]),
    (name="Benzene", formula="C6H6", category="Organic molecules", pg="D6h",
     geometry=benzene(1.397, 1.084), refs=[HERZBERG_1966]),
    (name="Methanol", formula="CH3OH", category="Organic molecules", pg="Cs",
     geometry=zmatrix([("C",), ("O", 1, 1.427), ("H", 2, 0.9561, 1, 108.87),
                       ("H", 1, 1.096, 2, 109.91, 3, 180.0), ("H", 1, 1.096, 2, 109.91, 3, 60.0),
                       ("H", 1, 1.096, 2, 109.91, 3, -60.0)]),
     refs=[("P. Venkateswarlu, W. Gordy, Methyl alcohol. II. Molecular structure, J. Chem. Phys. 23, 1200 (1955)",
            "10.1063/1.1742240"),
           ("E. Herbst, J. K. Messer, F. C. De Lucia, A new analysis and additional measurements of the millimeter " *
            "and submillimeter spectrum of methanol, J. Mol. Spectrosc. 108, 42 (1984)", "10.1016/0022-2852(84)90285-6")]),
    (name="Ethanol", formula="C2H5OH", category="Organic molecules", pg="Cs",
     geometry=zmatrix([("C",), ("C", 1, 1.512), ("O", 2, 1.4311, 1, 107.80), ("H", 3, 0.971, 2, 105.40, 1, 180.0),
                       ("H", 1, 1.098, 2, 107.81, 3, 180.0), ("H", 1, 1.088, 2, 107.82, 3, -58.90),
                       ("H", 1, 1.088, 2, 107.82, 3, 58.90), ("H", 2, 1.0859, 3, 107.83, 1, -121.10),
                       ("H", 2, 1.0859, 3, 107.83, 1, 121.10)]),
     refs=[("R. K. Kakar, C. R. Quade, Microwave rotational spectrum and internal rotation in gauche ethyl alcohol, " *
            "J. Chem. Phys. 72, 4300 (1980)", "10.1063/1.439723")],
     note="The trans (anti) conformer, a model geometry: check the parameters against the literature before " *
          "relying on them (the cited microwave study is of the gauche conformer)."),
    (name="Formaldehyde", formula="H2CO", category="Organic molecules", pg="C2v",
     geometry=zmatrix([("C",), ("O", 1, 1.205), ("H", 1, 1.111, 2, 121.93), ("H", 1, 1.111, 2, 121.93, 3, 180.0)]),
     refs=[HERZBERG_1966, GURVICH]),
    (name="Formic acid", formula="HCOOH", category="Organic molecules", pg="Cs",
     geometry=zmatrix([("C",), ("O", 1, 1.202), ("O", 1, 1.343, 2, 124.90), ("H", 1, 1.097, 3, 111.0, 2, 180.0),
                       ("H", 3, 0.972, 1, 106.3, 2, 0.0)]),
     refs=[HERZBERG_1966]),
    # Process gases and precursors (CVD, ALD, etching)
    (name="Silane", formula="SiH4", category="Process gases and precursors", pg="Td",
     geometry=tetrahedron("Si", "H", 1.4799),
     refs=[("D. R. J. Boyd, Infrared spectrum of trideuterosilane and the structure of the silane molecule, " *
            "J. Chem. Phys. 23, 922 (1955)", "10.1063/1.1742148")]),
    (name="Disilane", formula="Si2H6", category="Process gases and precursors", pg="D3d",
     geometry=staggered("Si", "H", 2.32, 1.47, 109.40),
     refs=[("G. W. Bethke, M. K. Wilson, Vibrational spectrum of disilane, J. Chem. Phys. 26, 1107 (1957)",
            "10.1063/1.1743480")],
     note="An approximate structure from the vibrational analysis (nearly tetrahedral angles)."),
    (name="Phosphine", formula="PH3", category="Process gases and precursors", pg="C3v",
     geometry=pyramid("P", "H", 1.421, 93.30), refs=[HERZBERG_1966]),
    (name="Arsine", formula="AsH3", category="Process gases and precursors", pg="C3v",
     geometry=pyramid("As", "H", 1.5108, 92.08),
     refs=[("G. S. Blevins, A. W. Jache, W. Gordy, Millimeter wave spectra of AsH3 and AsD3, Phys. Rev. 97, 684 " *
            "(1955)", "10.1103/PhysRev.97.684")]),
    (name="Diborane", formula="B2H6", category="Process gases and precursors", pg="D2h",
     geometry=diborane(1.76, 1.20, 1.317, 122.0), refs=[GURVICH]),
    (name="Dichlorosilane", formula="SiH2Cl2", category="Process gases and precursors", pg="C2v",
     geometry=xy2z2("Si", "H", "Cl", 1.480, 111.30, 2.033, 109.72),
     refs=[("R. W. Davis, M. C. L. Gerry, The microwave spectrum, structure, chlorine nuclear quadrupole coupling " *
            "constants, dipole moment, and centrifugal distortion constants of dichlorosilane, J. Mol. Spectrosc. " *
            "60, 117 (1976)", "10.1016/0022-2852(76)90120-X")]),
    (name="Silicon tetrachloride", formula="SiCl4", category="Process gases and precursors", pg="Td",
     geometry=tetrahedron("Si", "Cl", 2.0191),
     refs=[("R. R. Ryan, K. Hedberg, Effect of temperature on the structure of gaseous molecules. II. An " *
            "electron-diffraction investigation of the molecular structures of B2Cl4 and SiCl4, J. Chem. Phys. 50, " *
            "4986 (1969)", "10.1063/1.1670995")]),
    (name="Silicon tetrafluoride", formula="SiF4", category="Process gases and precursors", pg="Td",
     geometry=tetrahedron("Si", "F", 1.554),
     refs=[("M. Takami, H. Kuze, Infrared-microwave double resonance spectroscopy of the SiF4 ν3 fundamental " *
            "using a tunable diode laser, J. Chem. Phys. 78, 2204 (1983)", "10.1063/1.445063")]),
    (name="Nitrogen trifluoride", formula="NF3", category="Process gases and precursors", pg="C3v",
     geometry=pyramid("N", "F", 1.3648, 102.37),
     refs=[("S. E. Novick, W. Chen, M. R. Munrow, K. J. Grant, Hyperfine structure in the microwave spectrum of " *
            "NF3, J. Mol. Spectrosc. 179, 219 (1996)", "10.1006/jmsp.1996.0199"), KUCHITSU]),
    (name="Carbon tetrafluoride", formula="CF4", category="Process gases and precursors", pg="Td",
     geometry=tetrahedron("C", "F", 1.3151), refs=[KUCHITSU]),
    (name="Titanium tetrachloride", formula="TiCl4", category="Process gases and precursors", pg="Td",
     geometry=tetrahedron("Ti", "Cl", 2.1701),
     refs=[("Y. Morino, H. Uehara, Vibronic interactions in vanadium tetrachloride by gas electron diffraction, " *
            "J. Chem. Phys. 45, 4543 (1966)", "10.1063/1.1727535")],
     note="The cited electron-diffraction study is of VCl4; check that it reports this TiCl4 bond length."),
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
