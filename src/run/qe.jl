# Quantum ESPRESSO (pw.x) as a DFT engine: input from write_qe, results from
# the XML file pw.x writes at the end (energies, forces, stress, levels) and,
# for relaxations, every step from its output.

"""
Quantum ESPRESSO's pw.x. Found as `pw.x` on the PATH, or by the `command`
in the settings file (e.g. `"mpirun -np 4 pw.x"`). Pseudopotentials: the
*Pseudopotentials* field is a folder of UPF files, or `"PseudoDojo (as DFTK)"`:
the PseudoDojo files DFTK uses ([`dftk_pseudo_folder`](@ref)), so that the
two codes can be compared and the recommended cutoffs apply. The default is
`pseudo_dir` in the settings file, else `\$ESPRESSO_PSEUDO`, else
`~/espresso/pseudo` when it exists, else the PseudoDojo files. In a folder,
each element's file is the one whose name starts with the element symbol
followed by `.`, `_` or `-` (e.g. `Si.pbe-n-rrkjus_psl.1.0.0.UPF`,
`si_pbe_v1.uspp.F.UPF`).
"""
struct QEEngine <: DFTEngine end

"""The *Pseudopotentials* text that means DFTK's PseudoDojo files (for Quantum ESPRESSO)."""
const QE_DOJO = "PseudoDojo (as DFTK)"

const DFTK_PSEUDO_FOLDERS = Dict{String,String}()

"""
    dftk_pseudo_folder(family=DFTK_PSEUDO_FAMILY; env=default_dftk_env()) -> String

The folder of the UPF files of a PseudoPotentialData family, as DFTK
uses them (asked once per session from the DFTK environment, which
downloads the family the first time). Throws when DFTK is not installed.
"""
function dftk_pseudo_folder(family::AbstractString=DFTK_PSEUDO_FAMILY; env::AbstractString=default_dftk_env())
    haskey(DFTK_PSEUDO_FOLDERS, family) && return DFTK_PSEUDO_FOLDERS[family]
    dftk_env_ready(env) || throw(ArgumentError("the PseudoDojo files come with DFTK: install it (Run with DFTK once, " *
                                               "or install_dftk()), or give a folder of UPF files"))
    code = "using PseudoPotentialData; print(dirname(first(values(PseudoFamily(ARGS[1])))))"
    dir = strip(read(pipeline(julia_in(env, "-e", code, family); stderr=devnull), String))
    isdir(dir) || throw(ArgumentError("could not find the files of $family in $env"))
    return DFTK_PSEUDO_FOLDERS[family] = dir
end
engine_id(::QEEngine) = :qe
engine_name(::QEEngine) = "Quantum ESPRESSO"
engine_executables(::QEEngine) = ["pw.x"]
engine_install_hint(::QEEngine) = "e.g. conda install -c conda-forge qe"
engine_command(::QEEngine, command::Vector{String}, dir, name, input) = `$command -input $(basename(input))`
function engine_pseudo_default(e::QEEngine)
    s = engine_setting(e, "pseudo_dir")
    s isa AbstractString && !isempty(s) && return s
    haskey(ENV, "ESPRESSO_PSEUDO") && !isempty(ENV["ESPRESSO_PSEUDO"]) && return ENV["ESPRESSO_PSEUDO"]
    isdir(joinpath(homedir(), "espresso", "pseudo")) && return joinpath("~", "espresso", "pseudo")
    return QE_DOJO
end

"""
    qe_pseudo_files(dir, elements) -> Dict(element => file name)

For each element the UPF file in `dir` whose name starts with the element
symbol (any case) followed by `.`, `_` or `-`; when there are several, one
with the symbol's own case (`Si.…` before `si_…`), then the shortest. Throws an `ArgumentError` naming the elements without one.
"""
function qe_pseudo_files(dir::AbstractString, elements)
    isdir(dir) || throw(ArgumentError("no pseudopotential folder $dir (set it in the Pseudopotentials field)"))
    files = filter(f -> occursin(r"\.upf$"i, f), readdir(dir))
    out, missing_els = Dict{Symbol,String}(), Symbol[]
    for el in elements
        m = filter(f -> occursin(Regex("^" * string(el) * "[._-]", "i"), f), files)
        isempty(m) ? push!(missing_els, el) : (out[el] = sort(m; by=f -> (!startswith(f, string(el)), length(f)))[1])
    end
    isempty(missing_els) || throw(ArgumentError("no UPF file for " * join(missing_els, ", ") * " in $dir"))
    return out
end

function write_engine_inputs(::QEEngine, dir, name, c::Crystal; relax::Symbol, kspacing, ecutwfc, pseudo)
    pdir = strip(pseudo) == QE_DOJO ? dftk_pseudo_folder() : abspath(expanduser(pseudo))
    files = qe_pseudo_files(pdir, species_list(c))
    calc = relax === :none ? "scf" : relax === :atoms ? "relax" : "vc-relax"
    cc = relax === :atoms ? with_constraints(c; cell=(0, 0, 0, 0, 0, 0)) :
         relax === :cell && c.cell_constraints == (0, 0, 0, 0, 0, 0) ? with_constraints(c; cell=FREE_CELL) : c
    input = joinpath(dir, name * ".in")
    write_qe(input, cc; kspacing, ecutwfc, calculation=calc, prefix=name, pseudo_dir=pdir, outdir="./out",
             pseudo_files=files, run=true)
    return input, joinpath(dir, "out", name * ".xml")
end

# --- reading ------------------------------------------------------------------------

const RY_EV = HARTREE_EV / 2

"""The text of the first `<tag …>…</tag>` in `xml` (`nothing` if absent)."""
function xml_text(xml::AbstractString, tag::AbstractString)
    m = match(Regex("<" * tag * "(?:\\s[^>]*)?>(.*?)</" * tag * ">", "s"), xml)
    return m === nothing ? nothing : m[1]
end
xml_numbers(s) = s === nothing ? Float64[] : parse.(Float64, split(s))

"""
    read_qe_results(xml; log="", title="", clamped=nothing, start=nothing) -> DFTResults

The results of a pw.x run from its XML file (`outdir/prefix.xml`): final
structure, total energy, forces, Fermi level and the levels (eigenvalues,
occupations and weights per k-point and spin), with the relaxation steps
from the pw.x output `log` when it has several (starting from `start`, the
structure as written to the input; the step count is then theirs, which
includes the final SCF pw.x runs after a vc-relax).
"""
function read_qe_results(path::AbstractString; log::AbstractString="", title::AbstractString="", clamped=nothing,
                         start::Union{Nothing,Crystal}=nothing)
    xml = read(path, String)
    out = something(xml_text(xml, "output"), xml)
    st = xml_text(out, "atomic_structure")
    st === nothing && throw(ArgumentError("no atomic_structure in $path (not a pw.x XML file?)"))
    cell = [BOHR * Vec3(xml_numbers(xml_text(st, a))) for a in ("a1", "a2", "a3")]
    atoms = collect(eachmatch(r"<atom\s+name=\"([^\"]+)\"[^>]*>([^<]+)</atom>", st))
    species = [castep_species(m[1]) for m in atoms]
    pos = [BOHR * Vec3(parse.(Float64, split(m[2]))) for m in atoms]
    c = crystal_from_cartesian(Lattice(cell...), species, pos; title)
    energy = HARTREE_EV * only(xml_numbers(xml_text(xml_text(out, "total_energy"), "etot")))
    f = xml_numbers(xml_text(out, "forces"))
    forces = [Vec3(HA_PER_BOHR_EV_PER_A .* f[3i-2:3i]) for i in 1:length(f) ÷ 3]
    bs = xml_text(out, "band_structure")
    lsda = bs !== nothing && strip(something(xml_text(bs, "lsda"), "false")) == "true"
    nspin = lsda ? 2 : 1
    fermi = bs === nothing ? nothing : xml_text(bs, "fermi_energy")
    fermi === nothing && bs !== nothing && (fermi = xml_text(bs, "highestOccupiedLevel"))
    fermi = fermi === nothing ? nothing : HARTREE_EV * parse(Float64, strip(fermi))
    eigs, occs, weights = Vector{Float64}[], Vector{Float64}[], Float64[]
    if bs !== nothing
        ks = [m.match for m in eachmatch(r"<ks_energies>.*?</ks_energies>"s, bs)]
        w = [parse(Float64, match(r"weight=\"([^\"]+)\"", k)[1]) for k in ks]
        w ./= sum(w)
        e = [HARTREE_EV .* xml_numbers(xml_text(k, "eigenvalues")) for k in ks]
        o = [xml_numbers(xml_text(k, "occupations")) for k in ks]
        for σ in 1:nspin                                  # spin up for every k-point first, as DFTK
            for (i, k) in enumerate(ks)
                n = length(e[i]) ÷ nspin
                r = (σ - 1) * n + 1:σ * n
                push!(eigs, e[i][r]); push!(occs, (2 / nspin) .* o[i][r]); push!(weights, w[i])
            end
        end
    end
    conv = xml_text(out, "convergence_info")
    scf_ok = conv === nothing || strip(something(xml_text(xml_text(conv, "scf_conv"), "convergence_achieved"), "true")) == "true"
    opt = conv === nothing ? nothing : xml_text(conv, "opt_conv")
    opt_ok = opt === nothing || strip(something(xml_text(opt, "convergence_achieved"), "true")) == "true"
    nscf = conv === nothing ? 0 : parse(Int, strip(something(xml_text(xml_text(conv, "scf_conv"), "n_scf_steps"), "0")))
    nopt = opt === nothing ? 0 : parse(Int, strip(something(xml_text(opt, "n_opt_steps"), "0")))
    creator = match(r"<creator\s+NAME=\"[^\"]*\"\s+VERSION=\"([^\"]+)\"", xml)
    # the steps start from the structure the run started with (`start`, as written to the input), not the final one
    history = isempty(log) || !isfile(log) ? nothing : read_qe_history(log, something(start, c); clamped)
    history !== nothing && length(history) < 2 && (history = nothing)
    return DFTResults(; calculation=opt === nothing ? "scf" : "relax", converged=scf_ok && opt_ok,
                      iterations=opt === nothing ? nscf : history === nothing ? nopt : length(history) - 1, energy, structure=c, fermi, forces,
                      eigenvalues=eigs, occupations=occs, kweights=weights, n_spin=nspin,
                      version=creator === nothing ? "" : creator[1], history, program="Quantum ESPRESSO")
end

"""
    read_qe_history(log, like; clamped=nothing) -> RelaxHistory or nothing

Every SCF of a pw.x relax or vc-relax output: the energy ('!' line), the
forces, the stress (when printed) and the structure it was computed for (the
starting one from `like`, then each ATOMIC_POSITIONS / CELL_PARAMETERS block
printed after a step). Species and their order come from `like`, which must
list the atoms in the order of the input (species grouped, as `write_qe`
writes them).
"""
function read_qe_history(path::AbstractString, like::Crystal; clamped=nothing)
    lines = readlines(path)
    n = natoms(like)
    L = like.lattice
    current = like
    alat = NaN
    steps = NamedTuple[]
    final_block = false
    k = 0
    while (k += 1) <= length(lines)
        l = lines[k]
        if (m = match(r"lattice parameter \(alat\)\s*=\s*([-\d.Ee+]+)", l)) !== nothing
            alat = parse(Float64, m[1])
        elseif occursin("Begin final coordinates", l)
            final_block = true
        elseif occursin("End final coordinates", l)
            final_block = false
        elseif startswith(strip(l), "!") && occursin("total energy", l)
            e = parse(Float64, match(r"=\s*([-\d.Ee+]+)\s*Ry", l)[1]) * RY_EV
            push!(steps, (; energy=e, structure=current, forces=Vec3[], stress=nothing))
        elseif occursin("Forces acting on atoms", l) && !isempty(steps)
            fs = Vec3[]
            j = k
            while length(fs) < n && (j += 1) <= length(lines)
                m = match(r"atom\s+\d+\s+type\s+\d+\s+force\s*=\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)", lines[j])
                m === nothing || push!(fs, Vec3(parse.(Float64, (m[1], m[2], m[3]))) .* (RY_EV / BOHR))
            end
            steps[end] = merge(steps[end], (; forces=fs))
        elseif occursin(r"total\s+stress\s+\(Ry/bohr\*\*3\)", l) && !isempty(steps)
            σ = permutedims(reduce(hcat, [parse.(Float64, split(lines[k+i])[1:3]) for i in 1:3]))
            # pw.x prints −σ (positive = compressed); AtomBond's σ is positive when tensile; Ry/bohr³ → GPa
            steps[end] = merge(steps[end], (; stress=-σ .* (RY_EV / BOHR^3 * EV_PER_A3_IN_GPA)))
        elseif startswith(strip(l), "CELL_PARAMETERS") && !final_block
            m = match(r"alat\s*=\s*([-\d.Ee+]+)", l)
            scale = occursin(r"angstrom"i, l) ? 1.0 : occursin(r"bohr"i, l) ? BOHR :
                    BOHR * (m === nothing ? alat : parse(Float64, m[1]))
            L = Lattice((scale * Vec3(parse.(Float64, split(lines[k+i])[1:3])) for i in 1:3)...)
        elseif startswith(strip(l), "ATOMIC_POSITIONS") && !final_block
            rows = [split(lines[k+i]) for i in 1:n]
            xyz = [Vec3(parse.(Float64, r[2:4])) for r in rows]
            frac = occursin(r"crystal"i, l) ? xyz :
                   [cart_to_frac(L, (occursin(r"angstrom"i, l) ? 1.0 : occursin(r"bohr"i, l) ? BOHR : BOHR * alat) * p)
                    for p in xyz]
            current = Crystal(L, like.species, frac; pbc=like.pbc, title=like.title * " (QE step $(length(steps)))")
        end
    end
    steps = filter(s -> length(s.forces) == n, steps)
    # pw.x prints the stress of every SCF with tstress; it is the cell's only when the cell moves (vc-relax prints CELL_PARAMETERS)
    cell_moves = any(l -> startswith(strip(l), "CELL_PARAMETERS"), lines)
    cell_moves || (steps = [merge(s, (; stress=nothing)) for s in steps])
    isempty(steps) && return nothing
    cl = clamped === nothing ? fill(false, n) : collect(Bool, clamped)
    free = findall(!, cl)
    maxf = [isempty(free) ? 0.0 : maximum(norm, s.forces[free]) for s in steps]
    stresses = all(s -> s.stress !== nothing, steps) ? [s.stress for s in steps] : Matrix{Float64}[]
    return RelaxHistory([s.energy for s in steps], maxf, [s.structure for s in steps], [s.forces for s in steps], cl,
                        stresses)
end

"""The structure of a job as the engine writes it: atoms grouped by species (and which are clamped)."""
function grouped_like(c::Crystal)
    order = grouped_order(c)
    title = isempty(c.title) ? formula(c) : c.title
    return Crystal(c.lattice, c.species[order], c.frac[order]; pbc=c.pbc, title), [all(c.fixed[i]) for i in order]
end

function read_engine_results(::QEEngine, job)
    like, clamped = grouped_like(job.crystal)
    return read_qe_results(job.results; log=job.log, title=like.title, clamped, start=like)
end

function engine_history(::QEEngine, job)
    job.relax && isfile(job.log) || return nothing
    like, clamped = grouped_like(job.crystal)
    return try
        read_qe_history(job.log, like; clamped)
    catch                                          # a step half written: next time
        nothing
    end
end

register_engine!(QEEngine())
