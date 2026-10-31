"""Build AtomBond's structure library, data/library/, from its sources.

    python3 scripts/build_library.py [download_dir]

Crystals: CIF files from the Crystallography Open Database (COD,
https://www.crystallography.net), each a published experimental structure;
the CIF is stored as downloaded. III-V semiconductors: zincblende and wurtzite
cells written from the 300 K lattice constants of Vurgaftman, Meyer and
Ram-Mohan, J. Appl. Phys. 89, 5815 (2001) (one consistent review; COD's III-V
entries are mostly 1963 or older), with wurtzite u from Rietveld refinements. Molecules: experimental geometries from the
NIST Computational Chemistry Comparison and Benchmark Database (CCCBDB,
https://cccbdb.nist.gov, NIST Standard Reference Database 101), stored as XYZ.
Writes data/library/crystals/*.cif, data/library/molecules/*.xyz and the
index data/library/library.toml. To add a structure, add an entry to CRYSTALS
or MOLECULES and run the script again. Only the Python standard library is used.
"""

import html
import os
import re
import sys
import time
import urllib.request

# name, file, COD id, category, description, note
CRYSTALS = [
    dict(name="Silicon", file="Si.cif", cod=9011998, category="Group IV",
         description="Diamond structure (Fd-3m), a = 5.430941 Å at 25 °C; precision lattice constant"),
    dict(name="Diamond", file="C_diamond.cif", cod=9011997, category="Group IV",
         description="Diamond structure (Fd-3m), a = 3.566986 Å at 25 °C; precision lattice constant"),
    dict(name="Silicon-germanium Si0.5Ge0.5", file="SiGe.cif", cod=1537802, category="Group IV", ordered_sige=True,
         formula="SiGe",
         description="Ordered model of the random Si0.5Ge0.5 alloy: Si and Ge on the two diamond sublattices "
                     "(zincblende type, F-43m), measured alloy lattice constant a = 5.5373 Å",
         note="The real alloy is random (COD 1537802 has Si and Ge sharing one site 50/50). For DFT use this "
              "ordered cell, or a supercell with Si/Ge distributed as needed (e.g. an SQS)."),
    # Oxides
    dict(name="Quartz (alpha-SiO2)", file="SiO2_alpha_quartz.cif", cod=9013321, category="Oxides", formula="SiO2",
         description="Low (alpha) quartz, P3_221, a = 4.913437 Å, c = 5.405118 Å; high-resolution powder XRD"),
    dict(name="Corundum (alpha-Al2O3)", file="Al2O3_corundum.cif", cod=2300448, category="Oxides", formula="Al2O3",
         description="Corundum, R-3c (hexagonal setting), a = 4.758152 Å, c = 12.9897 Å at 298 K"),
    dict(name="Hafnia (monoclinic HfO2)", file="HfO2_monoclinic.cif", cod=9013470, category="Oxides", formula="HfO2",
         description="Monoclinic hafnia (baddeleyite type), P2_1/c, a = 5.1156, b = 5.1722, c = 5.2948 Å, beta = 99.18°"),
    dict(name="Zirconia (monoclinic ZrO2)", file="ZrO2_monoclinic.cif", cod=9016714, category="Oxides", formula="ZrO2",
         description="Baddeleyite, P2_1/c, a = 5.1505, b = 5.2116, c = 5.3173 Å, beta = 99.23°; neutron powder diffraction"),
    dict(name="Rutile (TiO2)", file="TiO2_rutile.cif", cod=9015662, category="Oxides", formula="TiO2",
         description="Rutile, P4_2/mnm, a = 4.5937 Å, c = 2.9587 Å; neutron powder diffraction"),
    dict(name="Anatase (TiO2)", file="TiO2_anatase.cif", cod=9015929, category="Oxides", formula="TiO2",
         description="Anatase, I4_1/amd, a = 3.7845 Å, c = 9.5143 Å; neutron powder diffraction"),
    dict(name="Zinc oxide (wurtzite ZnO)", file="ZnO_wurtzite.cif", cod=2300450, category="Oxides", formula="ZnO",
         description="Zincite, P6_3mc, a = 3.249308 Å, c = 5.205709 Å at 298 K"),
    dict(name="Magnesium oxide (MgO)", file="MgO_rocksalt.cif", cod=9005513, category="Oxides", formula="MgO",
         description="Periclase, rocksalt Fm-3m, a = 4.2113 Å"),
    dict(name="Gallium oxide (beta-Ga2O3)", file="Ga2O3_beta.cif", cod=2004987, category="Oxides", formula="Ga2O3",
         description="beta-Ga2O3, C2/m, a = 12.214, b = 3.0371, c = 5.7981 Å, beta = 103.83°"),
    dict(name="Strontium titanate (SrTiO3)", file="SrTiO3_perovskite.cif", cod=9006864, category="Oxides",
         formula="SrTiO3", description="Cubic perovskite, Pm-3m, a = 3.90528 Å (room temperature)"),
]

# III-V semiconductors: 300 K lattice constants from one review, so they are
# consistent with each other; wurtzite u from Rietveld refinements.
VURGAFTMAN = ("I. Vurgaftman, J. R. Meyer, L. R. Ram-Mohan, \"Band parameters for III-V compound semiconductors "
              "and their alloys\", J. Appl. Phys. 89, 5815 (2001)", "10.1063/1.1368156")
PASZKOWICZ_ALGAN = ("W. Paszkowicz, S. Podsiadlo, R. Minikayev, \"Rietveld-refinement study of aluminium and gallium "
                    "nitrides\", J. Alloys Compd. 382, 100 (2004)", "10.1016/j.jallcom.2004.05.036")
PASZKOWICZ_INN = ("W. Paszkowicz, R. Cerny, S. Krukowski, \"Rietveld refinement for indium nitride in the 105-295 K "
                  "range\", Powder Diffraction 18, 114 (2003)", "10.1154/1.1566957")
ZINCBLENDE = [("GaAs", "Gallium arsenide", 5.65325, "Table I"), ("AlAs", "Aluminium arsenide", 5.6611, "Table II"),
              ("InAs", "Indium arsenide", 6.0583, "Table III"), ("GaP", "Gallium phosphide", 5.4505, "Table IV"),
              ("AlP", "Aluminium phosphide", 5.4672, "Table V"), ("InP", "Indium phosphide", 5.8697, "Table VI"),
              ("GaSb", "Gallium antimonide", 6.0959, "Table VII"), ("AlSb", "Aluminium antimonide", 6.1355, "Table VIII"),
              ("InSb", "Indium antimonide", 6.4794, "Table IX")]
WURTZITE = [("GaN", "Gallium nitride", 3.189, 5.185, 0.3789, PASZKOWICZ_ALGAN,
             "Bulk values. Cross-check: NIST Structural Ceramics Database Z00623 (Leszczynski et al., "
             "J. Appl. Phys. 76, 4909 (1994)) gives a = 3.1880 Å, c = 5.18561 Å for a bulk single crystal at 294 K."),
            ("AlN", "Aluminium nitride", 3.112, 4.982, 0.3869, PASZKOWICZ_ALGAN, ""),
            ("InN", "Indium nitride", 3.545, 5.703, 0.3769, PASZKOWICZ_INN, "")]


def split_formula(f):
    return re.findall(r"[A-Z][a-z]?", f)


def prototype_cif(name, sg_symbol, sg_number, cell, sites, source):
    a, b, c, al, be, ga = cell
    lines = [f"# {name}, written by scripts/build_library.py", f"# {source}", f"data_{name.replace(' ', '_')}",
             f"_space_group_IT_number           {sg_number}", f"_symmetry_space_group_name_H-M   '{sg_symbol}'",
             f"_cell_length_a                   {a}", f"_cell_length_b                   {b}",
             f"_cell_length_c                   {c}", f"_cell_angle_alpha                {al}",
             f"_cell_angle_beta                 {be}", f"_cell_angle_gamma                {ga}",
             "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
             "_atom_site_fract_z"]
    lines += [f"{el}1 {el} {x} {y} {z}" for el, x, y, z in sites]
    return "\n".join(lines) + "\n"


def iii_v_entries():
    out = []
    for f, name, a, table in ZINCBLENDE:
        cat, an = split_formula(f)
        out.append(dict(name=f"{name} ({f})", file=f"{f}_zincblende.cif", formula=f, category="III-V",
                        cif=prototype_cif(f, "F -4 3 m", 216, (a, a, a, 90, 90, 90),
                                          [(cat, 0, 0, 0), (an, 0.25, 0.25, 0.25)], VURGAFTMAN[0]),
                        description=f"Zincblende, F-43m, a = {a} Å at 300 K",
                        source=f"Vurgaftman et al. 2001, {table} (recommended value)",
                        url="https://doi.org/" + VURGAFTMAN[1], reference=VURGAFTMAN[0], doi=VURGAFTMAN[1], note=""))
    for f, name, a, c, u, (uref, udoi), note in WURTZITE:
        cat, an = split_formula(f)
        out.append(dict(name=f"{name} ({f}, wurtzite)", file=f"{f}_wurtzite.cif", formula=f, category="III-V",
                        cif=prototype_cif(f, "P 63 m c", 186, (a, a, c, 90, 90, 120),
                                          [(cat, 0.3333333333, 0.6666666667, 0), (an, 0.3333333333, 0.6666666667, u)],
                                          f"{VURGAFTMAN[0]}; u: {uref}"),
                        description=f"Wurtzite, P6_3mc, a = {a} Å, c = {c} Å at 300 K, u = {u}",
                        source="Vurgaftman et al. 2001, Table X (a, c); Rietveld refinement (u)",
                        url="https://doi.org/" + VURGAFTMAN[1], reference=f"{VURGAFTMAN[0]}; u: {uref}",
                        doi=f"{VURGAFTMAN[1]}, {udoi}", note=note))
    return out


# name, file, CAS number (CCCBDB key), category, point group
# name, CAS number (the CCCBDB key, no dashes), formula, category, optional note.
# The XYZ file is named after the formula; the point group is read from CCCBDB.
SPIN_NOTE = "The ground state is a spin {0} (open shell): use a spin-polarized calculation."
MOLECULES = [
    # Diatomic molecules
    dict(name="Hydrogen", cas="1333740", formula="H2", category="Diatomic molecules"),
    dict(name="Nitrogen", cas="7727379", formula="N2", category="Diatomic molecules"),
    dict(name="Oxygen", cas="7782447", formula="O2", category="Diatomic molecules", note=SPIN_NOTE.format("triplet")),
    dict(name="Carbon monoxide", cas="630080", formula="CO", category="Diatomic molecules"),
    dict(name="Nitric oxide", cas="10102439", formula="NO", category="Diatomic molecules",
         note=SPIN_NOTE.format("doublet")),
    dict(name="Hydrogen fluoride", cas="7664393", formula="HF", category="Diatomic molecules"),
    dict(name="Hydrogen chloride", cas="7647010", formula="HCl", category="Diatomic molecules"),
    # Small molecules
    dict(name="Water", cas="7732185", formula="H2O", category="Small molecules"),
    dict(name="Ammonia", cas="7664417", formula="NH3", category="Small molecules"),
    dict(name="Carbon dioxide", cas="124389", formula="CO2", category="Small molecules"),
    dict(name="Hydrogen sulfide", cas="7783064", formula="H2S", category="Small molecules"),
    dict(name="Nitrous oxide", cas="10024972", formula="N2O", category="Small molecules"),
    dict(name="Sulfur dioxide", cas="7446095", formula="SO2", category="Small molecules"),
    dict(name="Ozone", cas="10028156", formula="O3", category="Small molecules"),
    dict(name="Hydrogen peroxide", cas="7722841", formula="H2O2", category="Small molecules"),
    dict(name="Hydrogen cyanide", cas="74908", formula="HCN", category="Small molecules"),
    # Organic molecules
    dict(name="Methane", cas="74828", formula="CH4", category="Organic molecules"),
    dict(name="Ethane", cas="74840", formula="C2H6", category="Organic molecules", point_group="D3d",
         note="CCCBDB labels this geometry D3h, but its coordinates are staggered ethane (D3d): they have an "
              "inversion centre and H-C-C-H dihedrals of 60 and 180 degrees."),
    dict(name="Ethylene", cas="74851", formula="C2H4", category="Organic molecules"),
    dict(name="Acetylene", cas="74862", formula="C2H2", category="Organic molecules"),
    dict(name="Benzene", cas="71432", formula="C6H6", category="Organic molecules"),
    dict(name="Methanol", cas="67561", formula="CH3OH", category="Organic molecules"),
    dict(name="Ethanol", cas="64175", formula="C2H5OH", category="Organic molecules"),
    dict(name="Formaldehyde", cas="50000", formula="H2CO", category="Organic molecules"),
    dict(name="Formic acid", cas="64186", formula="HCOOH", category="Organic molecules"),
    # Process gases and precursors (CVD, ALD, etching). CCCBDB has no experimental
    # geometry for GeH4, WF6, Ga(CH3)3 or Al(CH3)3 (checked October 2026).
    dict(name="Silane", cas="7803625", formula="SiH4", category="Process gases and precursors"),
    dict(name="Disilane", cas="1590870", formula="Si2H6", category="Process gases and precursors"),
    dict(name="Phosphine", cas="7803512", formula="PH3", category="Process gases and precursors"),
    dict(name="Arsine", cas="7784421", formula="AsH3", category="Process gases and precursors"),
    dict(name="Diborane", cas="19287457", formula="B2H6", category="Process gases and precursors"),
    dict(name="Dichlorosilane", cas="4109960", formula="SiH2Cl2", category="Process gases and precursors"),
    dict(name="Silicon tetrachloride", cas="10026047", formula="SiCl4", category="Process gases and precursors"),
    dict(name="Silicon tetrafluoride", cas="7783611", formula="SiF4", category="Process gases and precursors"),
    dict(name="Nitrogen trifluoride", cas="7783542", formula="NF3", category="Process gases and precursors"),
    dict(name="Carbon tetrafluoride", cas="75730", formula="CF4", category="Process gases and precursors"),
    dict(name="Titanium tetrachloride", cas="7550450", formula="TiCl4", category="Process gases and precursors"),
]


def get(url, path):
    if not os.path.exists(path):
        # A browser-like User-Agent (CCCBDB refuses Python's default one); one
        # request per second, and retries because COD resets connections under load.
        request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (AtomBond library builder)"})
        for attempt in range(5):
            try:
                time.sleep(1)
                with urllib.request.urlopen(request, timeout=60) as r, open(path, "wb") as f:
                    f.write(r.read())
                break
            except OSError:
                if os.path.exists(path):
                    os.remove(path)
                if attempt == 4:
                    raise
                time.sleep(3)
    with open(path, encoding="utf-8", errors="replace") as f:
        return f.read()


def cif_field(cif, name):
    m = re.search(rf"^{name}[ \t]+(?:'([^']*)'|(\S+))", cif, re.M)
    if m:
        return m.group(1) or m.group(2)
    m = re.search(rf"^{name}\s*\n;\s*\n(.*?)\n;", cif, re.M | re.S)
    return " ".join(m.group(1).split()) if m else ""


def cif_reference(cif):
    authors = re.search(r"_publ_author_name\s*\n((?:\s*'[^']*'\s*\n)+)", cif)
    names = re.findall(r"'([^']*)'", authors.group(1)) if authors else []
    who = ", ".join(names[:3]) + (" et al." if len(names) > 3 else "")
    title = cif_field(cif, "_publ_section_title")
    journal = cif_field(cif, "_journal_name_full")
    vol, page, year = (cif_field(cif, k) for k in ("_journal_volume", "_journal_page_first", "_journal_year"))
    return f"{who}, \"{title}\", {journal} {vol}, {page} ({year})", cif_field(cif, "_journal_paper_doi")


def ordered_sige(a, source_ref):
    return f"""# Ordered model of the random Si0.5Ge0.5 alloy, written by scripts/build_library.py.
# Si on 4a, Ge on 4c of F-43m (zincblende arrangement of the diamond lattice).
# Lattice constant: COD 1537802, {source_ref}
data_SiGe_ordered
_chemical_formula_sum            'Ge Si'
_space_group_IT_number           216
_symmetry_space_group_name_H-M   'F -4 3 m'
_cell_length_a                   {a}
_cell_length_b                   {a}
_cell_length_c                   {a}
_cell_angle_alpha                90
_cell_angle_beta                 90
_cell_angle_gamma                90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0 0 0
Ge1 Ge 0.25 0.25 0.25
"""


def cccbdb(page):
    """Release, Cartesian atoms and literature references of a CCCBDB experimental geometry page."""
    page = re.sub(r"<script.*?</script>", "", page, flags=re.S)
    lines = [l.strip() for l in html.unescape(re.sub(r"<[^>]+>", "\n", page)).splitlines() if l.strip()]
    release = re.match(r"Release \d+ \([^)]*\)", next(l for l in lines if l.startswith("Release"))).group(0)
    pg = re.search(r"Point Group (.*?)</H2>", page)
    point_group = html.unescape(re.sub(r"<[^>]+>", "", pg.group(1))).strip() if pg else ""
    if "z (Å)" not in lines:
        raise ValueError("no experimental Cartesian coordinates on the CCCBDB page")
    k = lines.index("z (Å)") + 1
    atoms = []
    while k + 3 < len(lines) and re.fullmatch(r"[A-Z][a-z]?\d+", lines[k]):
        atoms.append((re.match(r"[A-Z][a-z]?", lines[k]).group(0), *map(float, lines[k + 1:k + 4])))
        k += 4
    if not atoms:
        raise ValueError("no experimental coordinates at CCCBDB")
    refs, doi = [], []
    k = lines.index("reference") + 2 if "reference" in lines else len(lines)   # after "squib  reference  DOI"
    while k + 1 < len(lines) and not lines[k].startswith("Got a better number"):
        refs.append(lines[k + 1])
        k += 2
        if k < len(lines) and lines[k].startswith("10."):
            doi.append(lines[k])
            k += 1
    return release, point_group, atoms, refs, doi


def composition(formula):
    """Element counts of a formula such as "Ga(CH3)3"."""
    while (m := re.search(r"\(([A-Za-z0-9]+)\)(\d*)", formula)):
        formula = formula[:m.start()] + m.group(1) * int(m.group(2) or 1) + formula[m.end():]
    counts = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[el] = counts.get(el, 0) + int(n or 1)
    return counts


def toml_string(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def main():
    cache = sys.argv[1] if len(sys.argv) > 1 else os.path.join("scratch", "library_sources")
    os.makedirs(cache, exist_ok=True)
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "library")
    for sub in ("crystals", "molecules"):
        os.makedirs(os.path.join(root, sub), exist_ok=True)
    entries = []

    for c in CRYSTALS:
        url = f"https://www.crystallography.net/cod/{c['cod']}.cif"
        cif = get(url, os.path.join(cache, f"{c['cod']}.cif"))
        ref, doi = cif_reference(cif)
        out = ordered_sige(cif_field(cif, "_cell_length_a"), ref) if c.get("ordered_sige") else cif
        with open(os.path.join(root, "crystals", c["file"]), "w", encoding="utf-8") as f:
            f.write(out)
        entries.append(dict(name=c["name"], kind="crystal", category=c["category"],
                            formula=c.get("formula") or cif_field(cif, "_chemical_formula_sum").replace(" ", ""),
                            file="crystals/" + c["file"], description=c["description"],
                            source=f"Crystallography Open Database, COD {c['cod']}",
                            url=f"https://www.crystallography.net/cod/{c['cod']}.html",
                            reference=ref, doi=doi, note=c.get("note", "")))

    for e in iii_v_entries():
        with open(os.path.join(root, "crystals", e["file"]), "w", encoding="utf-8") as f:
            f.write(e.pop("cif"))
        entries.append(dict(name=e["name"], kind="crystal", category=e["category"], formula=e["formula"],
                            file="crystals/" + e.pop("file"), description=e["description"], source=e["source"],
                            url=e["url"], reference=e["reference"], doi=e["doi"], note=e["note"]))

    skipped = []
    for m in MOLECULES:
        url = f"https://cccbdb.nist.gov/expgeom2x.asp?casno={m['cas']}&charge=0"
        try:
            release, pg, atoms, refs, doi = cccbdb(get(url, os.path.join(cache, f"cccbdb_{m['cas']}.html")))
            got = {}
            for a in atoms:
                got[a[0]] = got.get(a[0], 0) + 1
            if got != composition(m["formula"]):        # CCCBDB shows another molecule for unknown entries
                raise ValueError("CCCBDB has no experimental geometry (it returned a different molecule)")
        except (ValueError, OSError, StopIteration) as err:
            skipped.append(f"{m['name']} ({err})")
            continue
        pg = m.get("point_group", pg)                  # corrections of CCCBDB's label
        cas = f"{m['cas'][:-3]}-{m['cas'][-3:-1]}-{m['cas'][-1]}"
        file = re.sub(r"[()]", "", m["formula"]) + ".xyz"
        with open(os.path.join(root, "molecules", file), "w", encoding="utf-8") as f:
            f.write(f"{len(atoms)}\n{m['name']} {m['formula']}, experimental geometry, NIST CCCBDB {release}, CAS {cas}\n")
            for s, x, y, z in atoms:
                f.write(f"{s:<2} {x:10.4f} {y:10.4f} {z:10.4f}\n")
        entries.append(dict(name=m["name"], kind="molecule", category=m["category"], formula=m["formula"],
                            file="molecules/" + file,
                            description="Experimental geometry" + (f", point group {pg}" if pg else ""),
                            source=f"NIST CCCBDB ({release}), CAS {cas}",
                            url=f"https://cccbdb.nist.gov/expgeom2x.asp?casno={m['cas']}",
                            reference="; ".join(refs), doi=", ".join(doi), note=m.get("note", "")))
    for s in skipped:
        print("skipped:", s)

    order = ["Group IV", "III-V", "Oxides", "Diatomic molecules", "Small molecules", "Organic molecules",
             "Process gases and precursors"]
    entries.sort(key=lambda e: order.index(e["category"]) if e["category"] in order else len(order))
    with open(os.path.join(root, "library.toml"), "w", encoding="utf-8") as f:
        f.write("# AtomBond structure library (generated by scripts/build_library.py; add structures there).\n"
                "# One [[structure]] per entry; `file` is relative to this folder.\n")
        for e in entries:
            f.write("\n[[structure]]\n")
            for k, val in e.items():
                f.write(f"{k} = {toml_string(val)}\n")
    print(f"{len(entries)} structures written to {os.path.normpath(root)}")


if __name__ == "__main__":
    main()
