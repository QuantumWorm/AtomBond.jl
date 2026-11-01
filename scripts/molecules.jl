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
    (name="Fluorine", formula="F2", category="Diatomic molecules", pg="D∞h",
     geometry=linear(("F", "F"), (1.41426,)),
     refs=[("T. R. Gilson, M. R. Purnell, Improved rotational constants for fluorine 19F2 using coherent anti-Stokes " *
            "Raman spectroscopy of the vibrational Q branch, J. Raman Spectrosc. 12, 222 (1982)", "10.1002/jrs.1250120304")]),
    (name="Chlorine", formula="Cl2", category="Diatomic molecules", pg="D∞h",
     geometry=linear(("Cl", "Cl"), (1.9882,)),
     refs=[("H. G. M. Edwards, E. A. M. Good, D. A. Long, Pure rotational Raman spectra of the chlorine species, 35Cl2 " *
            "and 35Cl37Cl, J. Chem. Soc., Faraday Trans. 2 72, 927 (1976)", "10.1039/F29767200927")],
     note="r_e of 35Cl2."),
    (name="Bromine", formula="Br2", category="Diatomic molecules", pg="D∞h",
     geometry=linear(("Br", "Br"), (2.28102,)),
     refs=[("C. Focsa, H. Li, P. F. Bernath, Characterization of the ground state of Br2 by laser-induced fluorescence " *
            "Fourier transform spectroscopy of the B3Π0+u–X1Σ+g system, J. Mol. Spectrosc. 200, 104 (2000)",
            "10.1006/jmsp.1999.8039")]),
    (name="Hydrogen bromide", formula="HBr", category="Diatomic molecules", pg="C∞v",
     geometry=linear(("H", "Br"), (1.4146,)),
     refs=[("G. Jones, W. Gordy, Submillimeter-wave spectra of HCl and HBr, Phys. Rev. 136, A1229 (1964)",
            "10.1103/PhysRev.136.A1229")]),
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
    (name="Carbonyl sulfide", formula="OCS", category="Small molecules", pg="C∞v",
     geometry=linear(("O", "C", "S"), (1.155386, 1.562021)),
     refs=[("A. Foord, J. G. Smith, D. H. Whiffen, An anharmonic force field for carbonyl sulphide, Mol. Phys. 29, " *
            "1685 (1975)", "10.1080/00268977500101491")]),
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
    # Etch and dopant gases
    (name="Boron trifluoride", formula="BF3", category="Process gases and precursors", pg="D3h",
     geometry=planar_xy3("B", "F", 1.30731),
     refs=[("A. G. Maki, T. Masiello, T. A. Blake, The high-resolution infrared spectrum of 10BF3 from 400 to " *
            "4600 cm-1, J. Mol. Struct. 742, 3 (2005)", "10.1016/j.molstruc.2004.11.089")]),
    (name="Boron trichloride", formula="BCl3", category="Process gases and precursors", pg="D3h",
     geometry=planar_xy3("B", "Cl", 1.7421),
     refs=[("S. Konaka, Y. Murata, K. Kuchitsu, Y. Morino, The molecular structure and force constants of boron " *
            "trifluoride and boron trichloride, Bull. Chem. Soc. Jpn. 39, 1134 (1966)", "10.1246/bcsj.39.1134")],
     note="Electron-diffraction r_g (thermal average), not an equilibrium r_e."),
    (name="Boron tribromide", formula="BBr3", category="Process gases and precursors", pg="D3h",
     geometry=planar_xy3("B", "Br", 1.896),
     refs=[("K.-G. Martinsen, N. Vogt, H. V. Volden, V. S. Lyutsarev, J. Vogt, Molecular structure and force field " *
            "of boron tribromide as determined from combined analysis of gas electron diffraction and spectroscopic " *
            "data and supported by quantum-chemical density-functional calculations, J. Mol. Struct. 385, 159 (1996)",
            "10.1016/S0022-2860(96)09435-5")],
     note="Equilibrium estimate r_he from electron diffraction and spectroscopy (r_g = 1.900 Å)."),
    (name="Sulfur hexafluoride", formula="SF6", category="Process gases and precursors", pg="Oh",
     geometry=octahedron("S", "F", 1.5560),
     refs=[("V. Boudon, J. L. Doménech, D. Bermejo, H. Willner, High-resolution Raman spectroscopy of the ν1 region " *
            "and Raman-Raman double resonance spectroscopy of the 2ν1-ν1 band of 32SF6 and 34SF6. Determination of " *
            "the equilibrium bond length of sulfur hexafluoride, J. Mol. Spectrosc. 228, 392 (2004)",
            "10.1016/j.jms.2004.04.006")]),
    (name="Fluoromethane", formula="CH3F", category="Process gases and precursors", pg="C3v",
     geometry=xy3z("C", "H", "F", 1.087, 110.2, 1.383),
     refs=[("J. Demaison, J. Breidung, W. Thiel, D. Papoušek, The equilibrium structure of methyl fluoride, " *
            "Struct. Chem. 10, 129 (1999)", "10.1023/A:1022085314343")],
     note="Semi-experimental r_e: experimental rotational constants combined with CCSD(T) values by least squares."),
    (name="Hexafluoroethane", formula="C2F6", category="Process gases and precursors", pg="D3d",
     geometry=staggered("C", "F", 1.56, 1.32, 109.5),
     refs=[("D. A. Swick, I. L. Karle, Structure and internal motion of C2F6 and Si2Cl6 vapors, J. Chem. Phys. 23, " *
            "1499 (1955)", "10.1063/1.1742337")],
     note="Early electron-diffraction structure of low precision (C–C ± 0.03 Å, C–C–F ± 1.5°): relax before use."),
    (name="Chlorine trifluoride", formula="ClF3", category="Process gases and precursors", pg="C2v",
     geometry=zmatrix([("Cl",), ("F", 1, 1.598), ("F", 1, 1.698, 2, 87.483), ("F", 1, 1.698, 2, 87.483, 3, 180.0)]),
     refs=[("D. F. Smith, The microwave spectrum and structure of chlorine trifluoride, J. Chem. Phys. 21, 609 (1953)",
            "10.1063/1.1698976")],
     note="Planar T shape; effective ground-state structure from the moments of inertia."),
    (name="Xenon difluoride", formula="XeF2", category="Process gases and precursors", pg="D∞h",
     geometry=linear(("F", "Xe", "F"), (1.9791, 1.9791)),
     refs=[("N. J. Brassington, H. G. M. Edwards, D. A. Long, Pure rotational Raman spectrum of XeF2, J. Chem. " *
            "Soc., Faraday Trans. 2 74, 1208 (1978)", "10.1039/F29787401208")],
     note="Ground-state r_0, not an equilibrium r_e."),
    (name="Germane", formula="GeH4", category="Process gases and precursors", pg="Td",
     geometry=tetrahedron("Ge", "H", 1.51714),
     refs=[("C. Richard, V. Boudon, A. Rizopoulos, J. Vander Auwera, F. Kwabia Tchana, Line positions and intensities " *
            "for the ν2/ν4 bands of 5 isotopologues of germane near 11.5 µm, J. Quant. Spectrosc. Radiat. Transfer " *
            "260, 107474 (2021)", "10.1016/j.jqsrt.2020.107474")]),
    (name="Phosphorus trifluoride", formula="PF3", category="Process gases and precursors", pg="C3v",
     geometry=pyramid("P", "F", 1.56099, 97.567),
     refs=[("H. Najib, Experimental rovibrational constants and equilibrium structure of phosphorus trifluoride, " *
            "J. Mol. Spectrosc. 305, 17 (2014)", "10.1016/j.jms.2014.09.008")]),
    (name="Phosphorus pentafluoride", formula="PF5", category="Process gases and precursors", pg="D3h",
     geometry=bipyramid("P", "F", 1.580, 1.532),
     refs=[("H. Kurimura, S. Yamamoto, T. Egawa, K. Kuchitsu, Geometrical structure of PF5 as studied by gas " *
            "electron diffraction and spectroscopic data, J. Mol. Struct. 140, 79 (1986)", "10.1016/0022-2860(86)80149-1")],
     note="Electron-diffraction r_g; the authors estimate r_e = 1.576 Å (axial) and 1.529 Å (equatorial)."),
    (name="Arsenic trifluoride", formula="AsF3", category="Process gases and precursors", pg="C3v",
     geometry=pyramid("As", "F", 1.7041, 95.767),
     refs=[("J. G. Smith, The anharmonic force field of arsenic trifluoride from the millimetre wave spectra of " *
            "vibrationally excited states, Mol. Phys. 35, 461 (1978)", "10.1080/00268977800100341")]),
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
