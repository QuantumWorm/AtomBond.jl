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
# Molecules: experimental gas-phase geometries (bond lengths and angles) from
# the cited spectroscopy and electron-diffraction literature; the coordinates
# are computed from those values (scripts/molecules.jl) and stored as XYZ.
# Writes data/library/crystals/*.cif, data/library/molecules/*.xyz and the
# index data/library/library.toml. To add a structure, add an entry to CRYSTALS
# or MOLECULES and run the script again. Only the Julia standard library is used.

using Downloads, Printf, LinearAlgebra

include("molecules.jl")
include(joinpath(@__DIR__, "..", "src", "build", "nanotube_geometry.jl"))

# name, file, COD id, category, description; optional formula, note, ordered_sige
const CRYSTALS = [
    (name="Silicon", file="Si.cif", cod=9011998, category="Semiconductors: group IV",
     description="Diamond structure (Fd-3m), a = 5.430941 Å at 25 °C; precision lattice constant"),
    (name="Diamond", file="C_diamond.cif", cod=9011997, category="Semiconductors: group IV",
     description="Diamond structure (Fd-3m), a = 3.566986 Å at 25 °C; precision lattice constant"),
    (name="Silicon-germanium Si0.5Ge0.5", file="SiGe.cif", cod=1537802, category="Semiconductors: group IV", ordered_sige=true,
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
    (name="Graphite", file="C_graphite.cif", cod=9011577, category="Carbon: graphite and graphene", formula="C",
     description="Hexagonal (2H, AB-stacked) graphite, P6_3/mmc, a = 2.464 Å, c = 6.711 Å; neutron diffraction on " *
                 "natural crystals"),
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
        push!(out, (name="$name ($f)", file="$(f)_zincblende.cif", formula=f, category="Semiconductors: III-V",
                    cif=prototype_cif(f, "F -4 3 m", 216, (a, a, a, 90, 90, 90),
                                      [(cat, 0, 0, 0), (an, 0.25, 0.25, 0.25)], VURGAFTMAN[1]),
                    description="Zincblende, F-43m, a = $a Å at 300 K",
                    source="Vurgaftman et al. 2001, $table (recommended value)",
                    url="https://doi.org/" * VURGAFTMAN[2], reference=VURGAFTMAN[1], doi=VURGAFTMAN[2], note=""))
    end
    for (f, name, a, c, u, (uref, udoi), note) in WURTZITE
        cat, an = split_formula(f)
        push!(out, (name="$name ($f, wurtzite)", file="$(f)_wurtzite.cif", formula=f, category="Semiconductors: III-V",
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


"""Text of `url`, downloaded to `path` unless it is already there."""
function get_text(url, path)
    if !isfile(path)
        # A browser-like User-Agent (some servers refuse default ones); one
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

const CRYSTAL_GROUPS = ["Semiconductors: group IV", "Semiconductors: III-V", "Oxides",
                        "Carbon: graphite and graphene", "Carbon: nanotubes"]

# Molecules: the usual chemical classification (organic, organometallic,
# inorganic) with groups by functional group or compound type, so a molecule
# is found where a chemistry textbook would put it. `use` (process use) is
# extra information, shown and searchable in the library tool.
const MOLECULE_GROUPS = ["Organic: hydrocarbons", "Organic: halocarbons", "Organic: alcohols and ethers",
                         "Organic: carbonyl compounds", "Organic: nitrogen and sulfur compounds",
                         "Organic: heterocycles", "Organic: radicals",
                         "Organometallic: metal alkyls", "Organometallic: metal carbonyls",
                         "Inorganic: elements", "Inorganic: hydrides", "Inorganic: oxides and sulfides", "Inorganic: oxoacids",
                         "Inorganic: nitrides", "Inorganic: nonmetal halides", "Inorganic: metal halides", "Inorganic: oxyhalides",
                         "Inorganic: radicals", "Carbon: fullerenes"]

# --- carbon nanostructures from graphite's in-plane lattice constant ------------------

"""Extended XYZ text: lattice vectors as columns of `L`, Cartesian positions, pbc flags."""
function extxyz(title, L, syms, pos, pbc)
    lat = join((@sprintf("%.6f", x) for x in vec(L)), " ")
    lines = ["$(length(syms))", "Lattice=\"$lat\" Properties=species:S:1:pos:R:3 pbc=\"$(join((p ? "T" : "F" for p in pbc), " "))\" $title"]
    append!(lines, [@sprintf("%-2s %12.6f %12.6f %12.6f", s, p...) for (s, p) in zip(syms, pos)])
    return join(lines, "\n") * "\n"
end

"""Graphene: one layer, periodic in-plane, `vacuum` Å along c (not periodic)."""
function graphene_xyz(a, vacuum)
    L = [a a/2 0; 0 a*sqrt(3)/2 0; 0 0 vacuum]
    pos = [L * [1/3, 1/3, 0.5], L * [2/3, 2/3, 0.5]]                # a1, a2 at 60°
    return extxyz("graphene", L, ["C", "C"], pos, (true, true, false))
end

"""
Single-wall (n, m) carbon nanotube from src/build/nanotube_geometry.jl (the
code behind AtomBond's `nanotube`): one translational period along c
(periodic), with `vacuum` Å between periodic images in a and b (not periodic).
"""
function nanotube_xyz(n, m, a, vacuum)
    pos, _, R, lt = nanotube_geometry(n, m, a / sqrt(3))
    box = 2R + vacuum
    pos = [p + [box / 2, box / 2, 0] for p in pos]
    return extxyz("($n,$m) nanotube", [box 0 0; 0 box 0; 0 0 lt], fill("C", length(pos)), pos, (false, false, true)), R, lt
end

const NANOTUBES = [(5, 5, "armchair, metallic"), (10, 10, "armchair, metallic"), (9, 0, "zigzag, small gap"),
                   (10, 0, "zigzag, semiconducting"), (6, 5, "chiral, semiconducting")]

"""Graphene and nanotube entries built from graphite's in-plane lattice constant `a`."""
function carbon_entries(a, cod, ref, doi)
    src = "Built from the graphite lattice constant a = $a Å (COD $cod)"
    common = (url="https://www.crystallography.net/cod/$cod.html", reference=ref, doi=doi)
    out = [(; name="Graphene", file="C_graphene.xyz", formula="C", category="Carbon: graphite and graphene",
              text=graphene_xyz(a, 15.0), description="Single layer, periodic in-plane, C–C = $(round(a / sqrt(3); digits=4)) Å, 15 Å vacuum along c",
              source=src, common..., note="One layer of the graphite structure (flat, as in graphite).")]
    for (n, m, kind) in NANOTUBES
        text, R, lt = nanotube_xyz(n, m, a, 10.0)
        push!(out, (; name="Carbon nanotube ($n,$m)", file="CNT_$(n)_$(m).xyz", formula="C", category="Carbon: nanotubes",
                    text, description="Single-wall ($n,$m) nanotube, $kind, diameter $(round(2R; digits=2)) Å, " *
                                      "one period $(round(lt; digits=3)) Å along c (periodic), 10 Å vacuum around it",
                    source=src, common...,
                    note="A model, not a measured structure: the graphene sheet rolled into a cylinder. Relax before use " *
                         "(curvature changes the bonds slightly, more for thin tubes)."))
    end
    return out
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

    g = only(c for c in CRYSTALS if c.file == "C_graphite.cif")
    gcif = get_text("https://www.crystallography.net/cod/$(g.cod).cif", joinpath(cache, "$(g.cod).cif"))
    for e in carbon_entries(parse(Float64, replace(cif_field(gcif, "_cell_length_a"), r"\(.*" => "")), g.cod,
                            cif_reference(gcif)...)
        write(joinpath(root, "crystals", e.file), e.text)
        push!(entries, (name=e.name, kind="crystal", category=e.category, formula=e.formula,
                        file="crystals/" * e.file, description=e.description, source=e.source,
                        url=e.url, reference=e.reference, doi=e.doi, note=e.note))
    end

    for e in iii_v_entries()
        write(joinpath(root, "crystals", e.file), e.cif)
        push!(entries, (name=e.name, kind="crystal", category=e.category, formula=e.formula,
                        file="crystals/" * e.file, description=e.description, source=e.source,
                        url=e.url, reference=e.reference, doi=e.doi, note=e.note))
    end

    for m in MOLECULES
        syms, pos = m.geometry
        file = strip(replace(m.formula, r"[()]" => "_"), '_') * ".xyz"   # Al(CH3)3 -> Al_CH3_3.xyz
        open(joinpath(root, "molecules", file), "w") do f
            println(f, length(syms))
            println(f, "$(m.name) $(m.formula), experimental geometry ($(m.refs[1][1]))")
            for (s, p) in zip(syms, orient(pos))
                @printf(f, "%-2s %10.4f %10.4f %10.4f\n", s, p...)
            end
        end
        m.category in MOLECULE_GROUPS || error("$(m.name): unknown group $(m.category)")
        dois = filter(!isempty, last.(m.refs))
        push!(entries, (name=m.name, kind="molecule", category=m.category, formula=m.formula,
                        file="molecules/" * file, description="Experimental geometry, point group $(m.pg)",
                        source="Experimental gas-phase structure from the literature",
                        url=isempty(dois) ? "" : "https://doi.org/" * first(dois),
                        reference=join(first.(m.refs), "; "), doi=join(dois, ", "), note=get(m, :note, ""),
                        use=get(m, :use, "")))
    end

    # Crystals keep their order within a group; molecules are sorted by name.
    order = [CRYSTAL_GROUPS; MOLECULE_GROUPS]
    rank(e) = something(findfirst(==(e.category), order), length(order) + 1)
    sort!(entries; by=e -> (rank(e), e.kind == "molecule" ? lowercase(e.name) : ""), alg=MergeSort)
    open(joinpath(root, "library.toml"), "w") do f
        println(f, "# AtomBond structure library (generated by scripts/build_library.jl; add structures there).")
        println(f, "# One [[structure]] per entry; `file` is relative to this folder.")
        for e in entries
            print(f, "\n[[structure]]\n")
            for (k, val) in pairs(e)
                k === :use && isempty(val) && continue
                println(f, "$k = $(toml_string(val))")
            end
        end
    end
    println("$(length(entries)) structures written to $root")
end

abspath(PROGRAM_FILE) == @__FILE__() && main()
