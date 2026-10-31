"""Build data/elements.csv, AtomBond's element table, from its published sources.

    python3 scripts/build_element_table.py [download_dir]

Downloads (or reads from download_dir, if the files are already there):

- NIST "Atomic Weights and Isotopic Compositions" (standard atomic weights;
  these are the IUPAC/CIAAW values),
  https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
- ASE's transcriptions of Cordero et al. covalent radii, Alvarez van der Waals
  radii, Jmol colors and element names, https://gitlab.com/ase/ase (ase/data)

Bondi van der Waals radii (with Mantina et al. for the main-group elements
Bondi omits) are kept from AtomBond's previous table (BONDI below).
Only the Python standard library is used.
"""

import os
import re
import sys
import urllib.request

SOURCES = {
    "nist.txt": "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some",
    "ase___init__.py": "https://gitlab.com/ase/ase/-/raw/master/ase/data/__init__.py",
    "ase_vdw_alvarez.py": "https://gitlab.com/ase/ase/-/raw/master/ase/data/vdw_alvarez.py",
    "ase_colors.py": "https://gitlab.com/ase/ase/-/raw/master/ase/data/colors.py",
}
ZMAX = 103

# IUPAC conventional atomic weights for elements whose standard atomic weight
# is an interval (CIAAW, abridged table).
CONVENTIONAL = {"H": 1.008, "Li": 6.94, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "Mg": 24.305,
                "Si": 28.085, "S": 32.06, "Cl": 35.45, "Ar": 39.95, "Br": 79.904, "Tl": 204.38, "Pb": 207.2}

# No standard atomic weight at NIST: mass number of the longest-lived isotope.
LONGEST_LIVED = {"Am": 243, "Cm": 247, "Bk": 247, "Cf": 251, "Es": 252, "Fm": 257, "Md": 258, "No": 259, "Lr": 266}

# Bondi, J. Phys. Chem. 68, 441 (1964); Mantina et al., J. Phys. Chem. A 113, 5806 (2009).
BONDI = {
    "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "Ne": 1.54, "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Ni": 1.63, "Cu": 1.40, "Zn": 1.39, "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Se": 1.90,
    "Br": 1.85, "Kr": 2.02, "Rb": 3.03, "Sr": 2.49, "Pd": 1.63, "Ag": 1.72, "Cd": 1.58, "In": 1.93, "Sn": 2.17,
    "Sb": 2.06, "Te": 2.06, "I": 1.98, "Xe": 2.16, "Cs": 3.43, "Ba": 2.68, "Pt": 1.72, "Au": 1.66, "Hg": 1.55,
    "Tl": 1.96, "Pb": 2.02, "Bi": 2.07, "Po": 1.97, "At": 2.02, "Rn": 2.20, "Fr": 3.48, "Ra": 2.83, "U": 1.86,
}


def fetch(folder):
    os.makedirs(folder, exist_ok=True)
    text = {}
    for name, url in SOURCES.items():
        path = os.path.join(folder, name)
        if not os.path.exists(path):
            urllib.request.urlretrieve(url, path)
        with open(path, encoding="utf-8") as f:
            text[name] = f.read()
    return text


def ase_array(src, name):
    """Entries of `name = np.array([...])` in an ASE data file, one per line, index 0 = dummy."""
    body = re.search(rf"^{name}\s*=\s*np\.array\(\s*\[(.*?)^\s*\]", src, re.S | re.M).group(1)
    rows = []
    for line in body.splitlines():
        line = line.split("#")[0].strip().rstrip(",")
        if line:
            rows.append(line)
    return rows


def number(s):
    s = s.strip()
    return None if s in ("np.nan", "missing", "_R_DUMMY") else float(s)


def weights(nist):
    out = {}
    for block in nist.split("Atomic Number = ")[1:]:
        z = int(block.split()[0])
        sym = re.search(r"Atomic Symbol = (\w+)", block).group(1)
        w = re.search(r"Standard Atomic Weight = (.*)", block).group(1).strip()
        if z in out or sym in ("D", "T"):
            continue
        if sym in LONGEST_LIVED:
            out[z] = (float(LONGEST_LIVED[sym]), "mass number")
        elif not re.match(r"[\d\[]", w):
            continue
        elif sym in CONVENTIONAL:
            out[z] = (CONVENTIONAL[sym], "")
        elif m := re.fullmatch(r"\[(\d+)\]", w):          # no stable isotope: mass number
            out[z] = (float(m.group(1)), "mass number")
        elif m := re.fullmatch(r"\[([\d.]+),([\d.]+)\]", w):
            out[z] = (round((float(m.group(1)) + float(m.group(2))) / 2, 4), "interval midpoint")
        else:
            out[z] = (float(re.match(r"[\d.]+", w).group(0)), "")
    return out


def main():
    folder = sys.argv[1] if len(sys.argv) > 1 else os.path.join("scratch", "element_sources")
    t = fetch(folder)
    init = t["ase___init__.py"]
    symbols = re.findall(r"'(\w+)'", re.search(r"chemical_symbols = \[(.*?)\]", init, re.S).group(1))
    names = re.findall(r"'(\w+)'", re.search(r"atomic_names = \[(.*?)\]", init, re.S).group(1))
    cov = [number(x) for x in ase_array(init, "covalent_radii")]
    vdw = [number(x) for x in ase_array(t["ase_vdw_alvarez.py"], "vdw_radii")]
    jmol = [tuple(float(c) for c in x.strip("()").split(",")) for x in ase_array(t["ase_colors.py"], "jmol_colors")]
    mass = weights(t["nist.txt"])

    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    with open(os.path.join(root, "data", "elements.csv"), "w", encoding="utf-8") as f:
        f.write(HEADER)
        f.write("Z,symbol,name,mass,mass_note,covalent_radius,vdw_radius_bondi,vdw_radius_alvarez,color_jmol\n")
        for z in range(1, ZMAX + 1):
            s = symbols[z]
            name = names[z - 1] if len(names) == len(symbols) - 1 else names[z]
            m, note = mass[z]
            r = cov[z] if z < len(cov) else None
            a = vdw[z] if z < len(vdw) else None
            col = "".join(f"{round(255 * c):02X}" for c in jmol[z]) if z < len(jmol) else ""
            fmt = lambda x: "" if x is None else f"{x:g}"
            f.write(f"{z},{s},{name},{m:.10g},{note},{fmt(r)},{fmt(BONDI.get(s))},{fmt(a)},{col}\n")


HEADER = """\
# AtomBond element table (generated by scripts/build_element_table.py; edit that
# script, or override values in ~/.atombond/elements.csv, rather than this file).
#
# mass                standard atomic weight (u), NIST "Atomic Weights and Isotopic
#                     Compositions" = IUPAC/CIAAW; IUPAC conventional value where the
#                     standard weight is an interval; mass_note "mass number" = no
#                     stable isotope, mass number of a long-lived isotope (NIST).
# covalent_radius     Å, B. Cordero et al., Dalton Trans. 2008, 2832 (C sp3, Mn/Fe/Co
#                     low spin); empty above Cm.
# vdw_radius_bondi    Å, A. Bondi, J. Phys. Chem. 68, 441 (1964), with M. Mantina et al.,
#                     J. Phys. Chem. A 113, 5806 (2009) for main-group elements Bondi omits.
# vdw_radius_alvarez  Å, S. Alvarez, Dalton Trans. 42, 8617 (2013); empty from Fm.
# color_jmol          Jmol's CPK-style element colors (https://jmol.sourceforge.net/jscolors/),
#                     the de facto standard: H white, C grey, N blue, O red, S yellow ...
#
# Lines starting with # are comments. Empty values mean "not tabulated".
"""

if __name__ == "__main__":
    main()
