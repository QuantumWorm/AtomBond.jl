# Build AtomBond's structure library, data/library/, from its sources.
#
#     julia scripts/build_library.jl [download_dir]
#
# Crystals: CIF files from the Crystallography Open Database (COD,
# https://www.crystallography.net), each a published experimental structure;
# the CIF is stored as downloaded. III-V semiconductors: zincblende and wurtzite
# cells written from the 300 K lattice constants of Vurgaftman, Meyer and
# Ram-Mohan, J. Appl. Phys. 89, 5815 (2001) (one consistent review; COD's III-V
# entries are mostly 1963 or older), with wurtzite u from Rietveld refinements.
# Molecules: experimental geometries from the NIST Computational Chemistry
# Comparison and Benchmark Database (CCCBDB, https://cccbdb.nist.gov, NIST
# Standard Reference Database 101), stored as XYZ.
# Writes data/library/crystals/*.cif, data/library/molecules/*.xyz and the
# index data/library/library.toml. To add a structure, add an entry to CRYSTALS
# or MOLECULES and run the script again. Only the Julia standard library is used.

using Downloads, Printf

# name, file, COD id, category, description; optional formula, note, ordered_sige
const CRYSTALS = [
    (name="Silicon", file="Si.cif", cod=9011998, category="Group IV",
     description="Diamond structure (Fd-3m), a = 5.430941 Å at 25 °C; precision lattice constant"),
    (name="Diamond", file="C_diamond.cif", cod=9011997, category="Group IV",
     description="Diamond structure (Fd-3m), a = 3.566986 Å at 25 °C; precision lattice constant"),
    (name="Silicon-germanium Si0.5Ge0.5", file="SiGe.cif", cod=1537802, category="Group IV", ordered_sige=true,
     formula="SiGe",
     description="Ordered model of the random Si0.5Ge0.5 alloy: Si and Ge on the two diamond sublattices " *
                 "(zincblende type, F-43m), measured alloy lattice constant a = 5.5373 Å",
     note="The real alloy is random (COD 1537802 has Si and Ge sharing one site 50/50). For DFT use this " *
          "ordered cell, or a supercell with Si/Ge distributed as needed (e.g. an SQS)."),
    # Oxides
    (name="Quartz (alpha-SiO2)", file="SiO2_alpha_quartz.cif", cod=9013321, category="Oxides", formula="SiO2",
     description="Low (alpha) quartz, P3_221, a = 4.913437 Å, c = 5.405118 Å; high-resolution powder XRD"),
    (name="Corundum (alpha-Al2O3)", file="Al2O3_corundum.cif", cod=2300448, category="Oxides", formula="Al2O3",
     description="Corundum, R-3c (hexagonal setting), a = 4.758152 Å, c = 12.9897 Å at 298 K"),
    (name="Hafnia (monoclinic HfO2)", file="HfO2_monoclinic.cif", cod=9013470, category="Oxides", formula="HfO2",
     description="Monoclinic hafnia (baddeleyite type), P2_1/c, a = 5.1156, b = 5.1722, c = 5.2948 Å, beta = 99.18°"),
    (name="Zirconia (monoclinic ZrO2)", file="ZrO2_monoclinic.cif", cod=9016714, category="Oxides", formula="ZrO2",
     description="Baddeleyite, P2_1/c, a = 5.1505, b = 5.2116, c = 5.3173 Å, beta = 99.23°; neutron powder diffraction"),
    (name="Rutile (TiO2)", file="TiO2_rutile.cif", cod=9015662, category="Oxides", formula="TiO2",
     description="Rutile, P4_2/mnm, a = 4.5937 Å, c = 2.9587 Å; neutron powder diffraction"),
    (name="Anatase (TiO2)", file="TiO2_anatase.cif", cod=9015929, category="Oxides", formula="TiO2",
     description="Anatase, I4_1/amd, a = 3.7845 Å, c = 9.5143 Å; neutron powder diffraction"),
    (name="Zinc oxide (wurtzite ZnO)", file="ZnO_wurtzite.cif", cod=2300450, category="Oxides", formula="ZnO",
     description="Zincite, P6_3mc, a = 3.249308 Å, c = 5.205709 Å at 298 K"),
    (name="Magnesium oxide (MgO)", file="MgO_rocksalt.cif", cod=9005513, category="Oxides", formula="MgO",
     description="Periclase, rocksalt Fm-3m, a = 4.2113 Å"),
    (name="Gallium oxide (beta-Ga2O3)", file="Ga2O3_beta.cif", cod=2004987, category="Oxides", formula="Ga2O3",
     description="beta-Ga2O3, C2/m, a = 12.214, b = 3.0371, c = 5.7981 Å, beta = 103.83°"),
    (name="Strontium titanate (SrTiO3)", file="SrTiO3_perovskite.cif", cod=9006864, category="Oxides",
     formula="SrTiO3", description="Cubic perovskite, Pm-3m, a = 3.90528 Å (room temperature)"),
]

# III-V semiconductors: 300 K lattice constants from one review, so they are
# consistent with each other; wurtzite u from Rietveld refinements.
const VURGAFTMAN = ("I. Vurgaftman, J. R. Meyer, L. R. Ram-Mohan, \"Band parameters for III-V compound semiconductors " *
                    "and their alloys\", J. Appl. Phys. 89, 5815 (2001)", "10.1063/1.1368156")
const PASZKOWICZ_ALGAN = ("W. Paszkowicz, S. Podsiadlo, R. Minikayev, \"Rietveld-refinement study of aluminium and gallium " *
                          "nitrides\", J. Alloys Compd. 382, 100 (2004)", "10.1016/j.jallcom.2004.05.036")
const PASZKOWICZ_INN = ("W. Paszkowicz, R. Cerny, S. Krukowski, \"Rietveld refinement for indium nitride in the 105-295 K " *
                        "range\", Powder Diffraction 18, 114 (2003)", "10.1154/1.1566957")
const ZINCBLENDE = [("GaAs", "Gallium arsenide", 5.65325, "Table I"), ("AlAs", "Aluminium arsenide", 5.6611, "Table II"),
                    ("InAs", "Indium arsenide", 6.0583, "Table III"), ("GaP", "Gallium phosphide", 5.4505, "Table IV"),
                    ("AlP", "Aluminium phosphide", 5.4672, "Table V"), ("InP", "Indium phosphide", 5.8697, "Table VI"),
                    ("GaSb", "Gallium antimonide", 6.0959, "Table VII"),
                    ("AlSb", "Aluminium antimonide", 6.1355, "Table VIII"),
                    ("InSb", "Indium antimonide", 6.4794, "Table IX")]
const WURTZITE = [("GaN", "Gallium nitride", 3.189, 5.185, 0.3789, PASZKOWICZ_ALGAN,
                   "Bulk values. Cross-check: NIST Structural Ceramics Database Z00623 (Leszczynski et al., " *
                   "J. Appl. Phys. 76, 4909 (1994)) gives a = 3.1880 Å, c = 5.18561 Å for a bulk single crystal at 294 K."),
                  ("AlN", "Aluminium nitride", 3.112, 4.982, 0.3869, PASZKOWICZ_ALGAN, ""),
                  ("InN", "Indium nitride", 3.545, 5.703, 0.3769, PASZKOWICZ_INN, "")]

split_formula(f) = [m.match for m in eachmatch(r"[A-Z][a-z]?", f)]

function prototype_cif(name, sg_symbol, sg_number, cell, sites, source)
    a, b, c, al, be, ga = cell
    lines = ["# $name, written by scripts/build_library.jl", "# $source", "data_$(replace(name, ' ' => '_'))",
             "_space_group_IT_number           $sg_number", "_symmetry_space_group_name_H-M   '$sg_symbol'",
             "_cell_length_a                   $a", "_cell_length_b                   $b",
             "_cell_length_c                   $c", "_cell_angle_alpha                $al",
             "_cell_angle_beta                 $be", "_cell_angle_gamma                $ga",
             "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
             "_atom_site_fract_z"]
    append!(lines, ["$(el)1 $el $x $y $z" for (el, x, y, z) in sites])
    return join(lines, "\n") * "\n"
end

"""Library entries (and CIF text) of the III-V semiconductors."""
function iii_v_entries()
    out = []
    for (f, name, a, table) in ZINCBLENDE
        cat, an = split_formula(f)
        push!(out, (name="$name ($f)", file="$(f)_zincblende.cif", formula=f, category="III-V",
                    cif=prototype_cif(f, "F -4 3 m", 216, (a, a, a, 90, 90, 90),
                                      [(cat, 0, 0, 0), (an, 0.25, 0.25, 0.25)], VURGAFTMAN[1]),
                    description="Zincblende, F-43m, a = $a Å at 300 K",
                    source="Vurgaftman et al. 2001, $table (recommended value)",
                    url="https://doi.org/" * VURGAFTMAN[2], reference=VURGAFTMAN[1], doi=VURGAFTMAN[2], note=""))
    end
    for (f, name, a, c, u, (uref, udoi), note) in WURTZITE
        cat, an = split_formula(f)
        push!(out, (name="$name ($f, wurtzite)", file="$(f)_wurtzite.cif", formula=f, category="III-V",
                    cif=prototype_cif(f, "P 63 m c", 186, (a, a, c, 90, 90, 120),
                                      [(cat, 0.3333333333, 0.6666666667, 0), (an, 0.3333333333, 0.6666666667, u)],
                                      "$(VURGAFTMAN[1]); u: $uref"),
                    description="Wurtzite, P6_3mc, a = $a Å, c = $c Å at 300 K, u = $u",
                    source="Vurgaftman et al. 2001, Table X (a, c); Rietveld refinement (u)",
                    url="https://doi.org/" * VURGAFTMAN[2], reference="$(VURGAFTMAN[1]); u: $uref",
                    doi="$(VURGAFTMAN[2]), $udoi", note))
    end
    return out
end

# name, CAS number (the CCCBDB key, no dashes), formula, category; optional
# note and point_group (a correction of CCCBDB's label).
# The XYZ file is named after the formula; the point group is read from CCCBDB.
spin_note(state) = "The ground state is a spin $state (open shell): use a spin-polarized calculation."
const MOLECULES = [
    # Diatomic molecules
    (name="Hydrogen", cas="1333740", formula="H2", category="Diatomic molecules"),
    (name="Nitrogen", cas="7727379", formula="N2", category="Diatomic molecules"),
    (name="Oxygen", cas="7782447", formula="O2", category="Diatomic molecules", note=spin_note("triplet")),
    (name="Carbon monoxide", cas="630080", formula="CO", category="Diatomic molecules"),
    (name="Nitric oxide", cas="10102439", formula="NO", category="Diatomic molecules", note=spin_note("doublet")),
    (name="Hydrogen fluoride", cas="7664393", formula="HF", category="Diatomic molecules"),
    (name="Hydrogen chloride", cas="7647010", formula="HCl", category="Diatomic molecules"),
    # Small molecules
    (name="Water", cas="7732185", formula="H2O", category="Small molecules"),
    (name="Ammonia", cas="7664417", formula="NH3", category="Small molecules"),
    (name="Carbon dioxide", cas="124389", formula="CO2", category="Small molecules"),
    (name="Hydrogen sulfide", cas="7783064", formula="H2S", category="Small molecules"),
    (name="Nitrous oxide", cas="10024972", formula="N2O", category="Small molecules"),
    (name="Sulfur dioxide", cas="7446095", formula="SO2", category="Small molecules"),
    (name="Ozone", cas="10028156", formula="O3", category="Small molecules"),
    (name="Hydrogen peroxide", cas="7722841", formula="H2O2", category="Small molecules"),
    (name="Hydrogen cyanide", cas="74908", formula="HCN", category="Small molecules"),
    # Organic molecules
    (name="Methane", cas="74828", formula="CH4", category="Organic molecules"),
    (name="Ethane", cas="74840", formula="C2H6", category="Organic molecules", point_group="D3d",
     note="CCCBDB labels this geometry D3h, but its coordinates are staggered ethane (D3d): they have an " *
          "inversion centre and H-C-C-H dihedrals of 60 and 180 degrees."),
    (name="Ethylene", cas="74851", formula="C2H4", category="Organic molecules"),
    (name="Acetylene", cas="74862", formula="C2H2", category="Organic molecules"),
    (name="Benzene", cas="71432", formula="C6H6", category="Organic molecules"),
    (name="Methanol", cas="67561", formula="CH3OH", category="Organic molecules"),
    (name="Ethanol", cas="64175", formula="C2H5OH", category="Organic molecules"),
    (name="Formaldehyde", cas="50000", formula="H2CO", category="Organic molecules"),
    (name="Formic acid", cas="64186", formula="HCOOH", category="Organic molecules"),
    # Process gases and precursors (CVD, ALD, etching). CCCBDB has no experimental
    # geometry for GeH4, WF6, Ga(CH3)3 or Al(CH3)3 (checked October 2026).
    (name="Silane", cas="7803625", formula="SiH4", category="Process gases and precursors"),
    (name="Disilane", cas="1590870", formula="Si2H6", category="Process gases and precursors"),
    (name="Phosphine", cas="7803512", formula="PH3", category="Process gases and precursors"),
    (name="Arsine", cas="7784421", formula="AsH3", category="Process gases and precursors"),
    (name="Diborane", cas="19287457", formula="B2H6", category="Process gases and precursors"),
    (name="Dichlorosilane", cas="4109960", formula="SiH2Cl2", category="Process gases and precursors"),
    (name="Silicon tetrachloride", cas="10026047", formula="SiCl4", category="Process gases and precursors"),
    (name="Silicon tetrafluoride", cas="7783611", formula="SiF4", category="Process gases and precursors"),
    (name="Nitrogen trifluoride", cas="7783542", formula="NF3", category="Process gases and precursors"),
    (name="Carbon tetrafluoride", cas="75730", formula="CF4", category="Process gases and precursors"),
    (name="Titanium tetrachloride", cas="7550450", formula="TiCl4", category="Process gases and precursors"),
]

"""Text of `url`, downloaded to `path` unless it is already there."""
function get_text(url, path)
    if !isfile(path)
        # A browser-like User-Agent (CCCBDB refuses some default ones); one
        # request per second, and retries because COD resets connections under load.
        for attempt in 1:5
            try
                sleep(1)
                Downloads.download(url, path; headers=["User-Agent" => "Mozilla/5.0 (AtomBond library builder)"],
                                   timeout=60)
                break
            catch
                rm(path; force=true)
                attempt == 5 && rethrow()
                sleep(3)
            end
        end
    end
    s = read(path, String)
    return isvalid(s) ? s : String(map(c -> isvalid(c) ? c : '�', collect(s)))
end

function cif_field(cif, name)
    m = match(Regex("^$name[ \\t]+(?:'([^']*)'|(\\S+))", "m"), cif)
    m === nothing || return something(m[1], m[2])
    m = match(Regex("^$name\\s*\\n;\\s*\\n(.*?)\\n;", "ms"), cif)
    return m === nothing ? "" : join(split(m[1]), " ")
end

function cif_reference(cif)
    authors = match(r"_publ_author_name\s*\n((?:\s*'[^']*'\s*\n)+)", cif)
    names = authors === nothing ? String[] : [m[1] for m in eachmatch(r"'([^']*)'", authors[1])]
    who = join(first(names, 3), ", ") * (length(names) > 3 ? " et al." : "")
    title = cif_field(cif, "_publ_section_title")
    journal = cif_field(cif, "_journal_name_full")
    vol, page, year = (cif_field(cif, k) for k in ("_journal_volume", "_journal_page_first", "_journal_year"))
    return "$who, \"$title\", $journal $vol, $page ($year)", cif_field(cif, "_journal_paper_doi")
end

ordered_sige(a, source_ref) = """
    # Ordered model of the random Si0.5Ge0.5 alloy, written by scripts/build_library.jl.
    # Si on 4a, Ge on 4c of F-43m (zincblende arrangement of the diamond lattice).
    # Lattice constant: COD 1537802, $source_ref
    data_SiGe_ordered
    _chemical_formula_sum            'Ge Si'
    _space_group_IT_number           216
    _symmetry_space_group_name_H-M   'F -4 3 m'
    _cell_length_a                   $a
    _cell_length_b                   $a
    _cell_length_c                   $a
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

# The HTML entities on CCCBDB pages; an unknown one stops the build rather
# than being copied into the library.
const ENTITIES = Dict("amp" => "&", "lt" => "<", "gt" => ">", "quot" => "\"", "apos" => "'", "nbsp" => " ",
                      "Aring" => "Å", "infin" => "∞", "rarr" => "→", "deg" => "°", "plusmn" => "±")

function unescape_html(s)
    return replace(s, r"&(#x[0-9A-Fa-f]+|#\d+|\w+);" => function (e)
        k = e[2:end-1]
        startswith(k, "#x") && return string(Char(parse(Int, k[3:end]; base=16)))
        startswith(k, "#") && return string(Char(parse(Int, k[2:end])))
        haskey(ENTITIES, k) || error("unknown HTML entity $e: add it to ENTITIES")
        return ENTITIES[k]
    end)
end

"""Release, Cartesian atoms and literature references of a CCCBDB experimental geometry page."""
function cccbdb(page)
    page = replace(page, r"<script.*?</script>"s => "")
    text = unescape_html(replace(page, r"<[^>]+>" => "\n"))
    lines = filter(!isempty, strip.(split(text, r"\R")))
    k = findfirst(l -> startswith(l, "Release"), lines)
    k === nothing && error("no CCCBDB release on the page")
    release = match(r"^Release \d+ \([^)]*\)", lines[k]).match
    pg = match(r"Point Group (.*?)</H2>", page)
    point_group = pg === nothing ? "" : strip(unescape_html(replace(pg[1], r"<[^>]+>" => "")))
    k = findfirst(==("z (Å)"), lines)
    k === nothing && error("no experimental Cartesian coordinates on the CCCBDB page")
    k += 1
    atoms = []
    while k + 3 <= length(lines) && occursin(r"^[A-Z][a-z]?\d+$", lines[k])
        push!(atoms, (match(r"^[A-Z][a-z]?", lines[k]).match, parse.(Float64, lines[k+1:k+3])...))
        k += 4
    end
    isempty(atoms) && error("no experimental coordinates at CCCBDB")
    refs, doi = String[], String[]
    r = findfirst(==("reference"), lines)
    k = r === nothing ? length(lines) + 1 : r + 2                 # after "squib  reference  DOI"
    while k + 1 <= length(lines) && !startswith(lines[k], "Got a better number")
        push!(refs, lines[k+1])
        k += 2
        if k <= length(lines) && startswith(lines[k], "10.")
            push!(doi, lines[k])
            k += 1
        end
    end
    return release, point_group, atoms, refs, doi
end

"""Element counts of a formula such as "Ga(CH3)3"."""
function composition(formula)
    while (m = match(r"\(([A-Za-z0-9]+)\)(\d*)", formula)) !== nothing
        n = isempty(m[2]) ? 1 : parse(Int, m[2])
        formula = formula[1:m.offset-1] * m[1]^n * formula[m.offset+ncodeunits(m.match):end]
    end
    counts = Dict{String,Int}()
    for m in eachmatch(r"([A-Z][a-z]?)(\d*)", formula)
        counts[m[1]] = get(counts, m[1], 0) + (isempty(m[2]) ? 1 : parse(Int, m[2]))
    end
    return counts
end

toml_string(s) = "\"" * replace(s, "\\" => "\\\\", "\"" => "\\\"") * "\""

function main(args=ARGS)
    cache = isempty(args) ? joinpath("scratch", "library_sources") : args[1]
    mkpath(cache)
    root = normpath(joinpath(@__DIR__, "..", "data", "library"))
    mkpath.(joinpath.(root, ("crystals", "molecules")))
    entries = []

    for c in CRYSTALS
        url = "https://www.crystallography.net/cod/$(c.cod).cif"
        cif = get_text(url, joinpath(cache, "$(c.cod).cif"))
        ref, doi = cif_reference(cif)
        out = get(c, :ordered_sige, false) ? ordered_sige(cif_field(cif, "_cell_length_a"), ref) : cif
        write(joinpath(root, "crystals", c.file), out)
        push!(entries, (name=c.name, kind="crystal", category=c.category,
                        formula=get(c, :formula, replace(cif_field(cif, "_chemical_formula_sum"), " " => "")),
                        file="crystals/" * c.file, description=c.description,
                        source="Crystallography Open Database, COD $(c.cod)",
                        url="https://www.crystallography.net/cod/$(c.cod).html",
                        reference=ref, doi, note=get(c, :note, "")))
    end

    for e in iii_v_entries()
        write(joinpath(root, "crystals", e.file), e.cif)
        push!(entries, (name=e.name, kind="crystal", category=e.category, formula=e.formula,
                        file="crystals/" * e.file, description=e.description, source=e.source,
                        url=e.url, reference=e.reference, doi=e.doi, note=e.note))
    end

    skipped = String[]
    for m in MOLECULES
        url = "https://cccbdb.nist.gov/expgeom2x.asp?casno=$(m.cas)&charge=0"
        local release, pg, atoms, refs, doi
        try
            release, pg, atoms, refs, doi = cccbdb(get_text(url, joinpath(cache, "cccbdb_$(m.cas).html")))
            got = Dict{String,Int}()
            for a in atoms
                got[a[1]] = get(got, a[1], 0) + 1
            end
            got == composition(m.formula) ||          # CCCBDB shows another molecule for unknown entries
                error("CCCBDB has no experimental geometry (it returned a different molecule)")
        catch err
            push!(skipped, "$(m.name) ($(sprint(showerror, err)))")
            continue
        end
        pg = get(m, :point_group, pg)                  # corrections of CCCBDB's label
        cas = "$(m.cas[1:end-3])-$(m.cas[end-2:end-1])-$(m.cas[end])"
        file = replace(m.formula, r"[()]" => "") * ".xyz"
        open(joinpath(root, "molecules", file), "w") do f
            println(f, length(atoms))
            println(f, "$(m.name) $(m.formula), experimental geometry, NIST CCCBDB $release, CAS $cas")
            for (s, x, y, z) in atoms
                @printf(f, "%-2s %10.4f %10.4f %10.4f\n", s, x, y, z)
            end
        end
        push!(entries, (name=m.name, kind="molecule", category=m.category, formula=m.formula,
                        file="molecules/" * file,
                        description="Experimental geometry" * (isempty(pg) ? "" : ", point group $pg"),
                        source="NIST CCCBDB ($release), CAS $cas",
                        url="https://cccbdb.nist.gov/expgeom2x.asp?casno=$(m.cas)",
                        reference=join(refs, "; "), doi=join(doi, ", "), note=get(m, :note, "")))
    end
    foreach(s -> println("skipped: ", s), skipped)

    order = ["Group IV", "III-V", "Oxides", "Diatomic molecules", "Small molecules", "Organic molecules",
             "Process gases and precursors"]
    sort!(entries; by=e -> something(findfirst(==(e.category), order), length(order) + 1), alg=MergeSort)
    open(joinpath(root, "library.toml"), "w") do f
        println(f, "# AtomBond structure library (generated by scripts/build_library.jl; add structures there).")
        println(f, "# One [[structure]] per entry; `file` is relative to this folder.")
        for e in entries
            print(f, "\n[[structure]]\n")
            for (k, val) in pairs(e)
                println(f, "$k = $(toml_string(val))")
            end
        end
    end
    println("$(length(entries)) structures written to $root")
end

abspath(PROGRAM_FILE) == @__FILE__() && main()
