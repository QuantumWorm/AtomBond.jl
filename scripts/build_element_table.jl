# Build data/elements.csv, AtomBond's element table, from its published sources.
#
#     julia scripts/build_element_table.jl [download_dir]
#
# Downloads (or reads from download_dir, if the files are already there):
#
# - NIST "Atomic Weights and Isotopic Compositions" (standard atomic weights;
#   these are the IUPAC/CIAAW values),
#   https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
# - ASE's transcriptions of Cordero et al. covalent radii, Alvarez van der Waals
#   radii, Jmol colors and element names, https://gitlab.com/ase/ase (ase/data)
#
# Bondi van der Waals radii (with Mantina et al. for the main-group elements
# Bondi omits) are kept from AtomBond's previous table (BONDI below).
# Only the Julia standard library is used.

using Downloads, Printf

const SOURCES = [
    "nist.txt" => "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some",
    "ase___init__.py" => "https://gitlab.com/ase/ase/-/raw/master/ase/data/__init__.py",
    "ase_vdw_alvarez.py" => "https://gitlab.com/ase/ase/-/raw/master/ase/data/vdw_alvarez.py",
    "ase_colors.py" => "https://gitlab.com/ase/ase/-/raw/master/ase/data/colors.py",
]
const ZMAX = 103

# IUPAC conventional atomic weights for elements whose standard atomic weight
# is an interval (CIAAW, abridged table).
const CONVENTIONAL = Dict("H" => 1.008, "Li" => 6.94, "B" => 10.81, "C" => 12.011, "N" => 14.007, "O" => 15.999,
                          "Mg" => 24.305, "Si" => 28.085, "S" => 32.06, "Cl" => 35.45, "Ar" => 39.95, "Br" => 79.904,
                          "Tl" => 204.38, "Pb" => 207.2)

# No standard atomic weight at NIST: mass number of the longest-lived isotope.
const LONGEST_LIVED = Dict("Am" => 243, "Cm" => 247, "Bk" => 247, "Cf" => 251, "Es" => 252, "Fm" => 257,
                           "Md" => 258, "No" => 259, "Lr" => 266)

# Bondi, J. Phys. Chem. 68, 441 (1964); Mantina et al., J. Phys. Chem. A 113, 5806 (2009).
const BONDI = Dict(
    "H" => 1.20, "He" => 1.40, "Li" => 1.82, "Be" => 1.53, "B" => 1.92, "C" => 1.70, "N" => 1.55, "O" => 1.52,
    "F" => 1.47, "Ne" => 1.54, "Na" => 2.27, "Mg" => 1.73, "Al" => 1.84, "Si" => 2.10, "P" => 1.80, "S" => 1.80,
    "Cl" => 1.75, "Ar" => 1.88, "K" => 2.75, "Ca" => 2.31, "Ni" => 1.63, "Cu" => 1.40, "Zn" => 1.39, "Ga" => 1.87,
    "Ge" => 2.11, "As" => 1.85, "Se" => 1.90, "Br" => 1.85, "Kr" => 2.02, "Rb" => 3.03, "Sr" => 2.49, "Pd" => 1.63,
    "Ag" => 1.72, "Cd" => 1.58, "In" => 1.93, "Sn" => 2.17, "Sb" => 2.06, "Te" => 2.06, "I" => 1.98, "Xe" => 2.16,
    "Cs" => 3.43, "Ba" => 2.68, "Pt" => 1.72, "Au" => 1.66, "Hg" => 1.55, "Tl" => 1.96, "Pb" => 2.02, "Bi" => 2.07,
    "Po" => 1.97, "At" => 2.02, "Rn" => 2.20, "Fr" => 3.48, "Ra" => 2.83, "U" => 1.86,
)

function fetch(folder)
    mkpath(folder)
    return Dict(map(SOURCES) do (name, url)
        path = joinpath(folder, name)
        isfile(path) || Downloads.download(url, path)
        name => read(path, String)
    end)
end

"""Entries of `name = np.array([...])` in an ASE data file, one per line, index 1 = dummy."""
function ase_array(src, name)
    body = match(Regex("^$name\\s*=\\s*np\\.array\\(\\s*\\[(.*?)^\\s*\\]", "sm"), src)[1]
    rows = String[]
    for line in split(body, '\n')
        line = rstrip(strip(split(line, '#')[1]), ',')
        isempty(line) || push!(rows, line)
    end
    return rows
end

number(s) = (s = strip(s); s in ("np.nan", "missing", "_R_DUMMY") ? nothing : parse(Float64, s))

function weights(nist)
    out = Dict{Int,Tuple{Float64,String}}()
    for block in split(nist, "Atomic Number = ")[2:end]
        z = parse(Int, split(block)[1])
        sym = match(r"Atomic Symbol = (\w+)", block)[1]
        w = strip(match(r"Standard Atomic Weight = (.*)", block)[1])
        (haskey(out, z) || sym in ("D", "T")) && continue
        if haskey(LONGEST_LIVED, sym)
            out[z] = (float(LONGEST_LIVED[sym]), "mass number")
        elseif !occursin(r"^[\d\[]", w)
            continue
        elseif haskey(CONVENTIONAL, sym)
            out[z] = (CONVENTIONAL[sym], "")
        elseif (m = match(r"^\[(\d+)\]$", w)) !== nothing          # no stable isotope: mass number
            out[z] = (parse(Float64, m[1]), "mass number")
        elseif (m = match(r"^\[([\d.]+),([\d.]+)\]$", w)) !== nothing
            out[z] = (round((parse(Float64, m[1]) + parse(Float64, m[2])) / 2; digits=4), "interval midpoint")
        else
            out[z] = (parse(Float64, match(r"^[\d.]+", w).match), "")
        end
    end
    return out
end

quoted(list) = [m[1] for m in eachmatch(r"'(\w+)'", list)]

function main(args=ARGS)
    folder = isempty(args) ? joinpath("scratch", "element_sources") : args[1]
    t = fetch(folder)
    init = t["ase___init__.py"]
    symbols = quoted(match(r"chemical_symbols = \[(.*?)\]"s, init)[1])
    names = quoted(match(r"atomic_names = \[(.*?)\]"s, init)[1])
    cov = number.(ase_array(init, "covalent_radii"))
    vdw = number.(ase_array(t["ase_vdw_alvarez.py"], "vdw_radii"))
    jmol = [parse.(Float64, split(strip(x, ['(', ')']), ',')) for x in ase_array(t["ase_colors.py"], "jmol_colors")]
    mass = weights(t["nist.txt"])

    # The ASE lists start with a dummy element (index 1 here is Z = 0).
    fmt(x) = x === nothing ? "" : @sprintf("%g", x)
    open(joinpath(@__DIR__, "..", "data", "elements.csv"), "w") do f
        print(f, HEADER)
        println(f, "Z,symbol,name,mass,mass_note,covalent_radius,vdw_radius_bondi,vdw_radius_alvarez,color_jmol")
        for z in 1:ZMAX
            s = symbols[z+1]
            name = length(names) == length(symbols) - 1 ? names[z] : names[z+1]
            m, note = mass[z]
            r = z < length(cov) ? cov[z+1] : nothing
            a = z < length(vdw) ? vdw[z+1] : nothing
            col = z < length(jmol) ? join(@sprintf("%02X", round(Int, 255c)) for c in jmol[z+1]) : ""
            println(f, join((z, s, name, @sprintf("%.10g", m), note, fmt(r), fmt(get(BONDI, s, nothing)), fmt(a), col), ","))
        end
    end
end

const HEADER = """
# AtomBond element table (generated by scripts/build_element_table.jl; edit that
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
# Covalent radii, Alvarez radii and Jmol colours were transcribed from ASE's data
# files (https://gitlab.com/ase/ase, ase/data), which cite the same sources.
#
# Lines starting with # are comments. Empty values mean "not tabulated".
"""

abspath(PROGRAM_FILE) == @__FILE__() && main()
