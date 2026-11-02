# Running DFTK.jl from AtomBond: the calculation runs in a separate Julia
# process, in its own environment (DFTK is large and not a dependency of
# AtomBond), from the script that write_dftk writes. The window stays
# responsive; the job's log file shows the progress and a TOML file holds the
# results when it is done.

import TOML

"""Packages of the DFTK environment."""
const DFTK_PACKAGES = ["DFTK", "PseudoPotentialData", "AtomsBase", "Unitful", "UnitfulAtomic", "GeometryOptimization"]

"""
    default_dftk_env() -> String

The Julia environment AtomBond runs DFTK in: `environments/atombond-dftk`
in the first depot (e.g. `~/.julia/environments/atombond-dftk`), shared by
all AtomBond sessions. [`install_dftk`](@ref) creates it.
"""
default_dftk_env() = joinpath(first(DEPOT_PATH), "environments", "atombond-dftk")

"""Whether `env` is an environment with DFTK and the other packages the scripts need."""
function dftk_env_ready(env::AbstractString=default_dftk_env())
    project = joinpath(env, "Project.toml")
    isfile(project) || return false
    deps = try
        get(TOML.parsefile(project), "deps", Dict())
    catch
        return false
    end
    return all(p -> haskey(deps, p), DFTK_PACKAGES)
end

"""
A Julia command for a helper process: this Julia, no startup file, no colours
in its output (the log is shown as text), `env` as the project.
"""
function julia_in(env::AbstractString, args...; dir::AbstractString=pwd())
    cmd = `$(Base.julia_cmd()) --startup-file=no --color=no --project=$env $args`
    # Not the caller's load path or project (e.g. the test environment of Pkg.test).
    vars = filter(((k, _),) -> !(k in ("JULIA_LOAD_PATH", "JULIA_PROJECT")), ENV)
    return Cmd(setenv(cmd, vars); dir)
end

"""
    DFTKJob

A DFTK calculation started by [`run_dftk`](@ref) (or the DFTK installation
started by [`install_dftk`](@ref)): the process, its folder, script, log and
results files, and the structure it was started from.
"""
mutable struct DFTKJob
    crystal::Union{Nothing,Crystal}
    dir::String
    script::String
    log::String
    results::String
    process::Base.Process
    started::Float64
    stopped::Bool
end

"""
    install_dftk(env=default_dftk_env()) -> DFTKJob

Create the DFTK environment (adds DFTK, PseudoPotentialData, AtomsBase,
Unitful, UnitfulAtomic and GeometryOptimization, and precompiles them) in a
separate process. It downloads and compiles a lot, so it takes several
minutes the first time; [`job_status`](@ref) says when it is done and the log
shows the progress.
"""
function install_dftk(env::AbstractString=default_dftk_env())
    mkpath(env)
    log = joinpath(env, "install.log")
    code = "using Pkg; Pkg.add($(repr(DFTK_PACKAGES))); Pkg.precompile(); println(\"DFTK environment ready\")"
    proc = open(log, "w") do io
        run(pipeline(julia_in(env, "-e", code; dir=env); stdout=io, stderr=io); wait=false)
    end
    return DFTKJob(nothing, env, "", log, "", proc, time(), false)
end

"""
    run_dftk(c; dir, name, env=default_dftk_env(), kspacing=0.25, ecutwfc=40, pseudo) -> DFTKJob

Start a DFTK calculation of `c` in the background: writes `name.jl` (see
[`write_dftk`](@ref): an SCF, or a relaxation when `c` has constraints) into
`dir` and runs it with Julia in `env`, logging to `name.log`. When it
finishes, `name.results.toml` holds the results ([`job_results`](@ref)).
`dir` defaults to a new folder `dftk_<name>_<time>` under the current folder.
With `bands=true` an SCF is followed by a non-self-consistent step: the
levels on a k-grid twice as dense (a fuller DOS) and, for a crystal without
vacuum, the band structure along the high-symmetry path, in
`name.bands.toml` ([`read_dftk_results`](@ref) reads it too).
"""
function run_dftk(c::Crystal; name::AbstractString=formula(c),
                  dir::AbstractString=joinpath(pwd(), "dftk_" * name * "_" * Dates.format(Dates.now(), "yyyymmdd-HHMMSS")),
                  env::AbstractString=default_dftk_env(), kspacing::Real=0.25, ecutwfc::Real=40,
                  pseudo::AbstractString=DFTK_PSEUDO_FAMILY, bands::Bool=false)
    dftk_env_ready(env) || throw(ArgumentError("no DFTK environment at $env; create it with install_dftk()"))
    occursin(r"^[\w.+-]+$", name) || throw(ArgumentError("use letters, digits, . + - _ in the name"))
    mkpath(dir)
    script, log, results = (joinpath(dir, name * ext) for ext in (".jl", ".log", ".results.toml"))
    bandsfile = joinpath(dir, name * ".bands.toml")
    foreach(f -> isfile(f) && rm(f), (results, bandsfile))
    write_dftk(script, c; kspacing, ecutwfc, pseudo, results=basename(results),
               bands=bands ? basename(bandsfile) : nothing)
    proc = open(log, "w") do io
        run(pipeline(julia_in(env, script; dir); stdout=io, stderr=io); wait=false)
    end
    return DFTKJob(c, dir, script, log, results, proc, time(), false)
end

"""
    job_status(job) -> Symbol

`:running`, `:done` (finished and, for a calculation, wrote its results),
`:failed` (see the log) or `:stopped` ([`stop_job!`](@ref)).
"""
function job_status(job::DFTKJob)
    process_running(job.process) && return :running
    job.stopped && return :stopped
    ok = success(job.process) && (isempty(job.results) || isfile(job.results))
    return ok ? :done : :failed
end

# --- recommended cutoffs ------------------------------------------------------------

"""Cutoff hints (Ha) per (family, element), filled by [`dftk_cutoff_hints`](@ref)."""
const CUTOFF_HINTS = Dict{Tuple{String,Symbol},NamedTuple{(:low, :normal, :high),NTuple{3,Float64}}}()

"""
    dftk_cutoff_hints(elements; family=DFTK_PSEUDO_FAMILY, env=default_dftk_env()) -> Dict

The kinetic-energy cutoffs (Ha) the pseudopotential family recommends for
each element, `element => (; low, normal, high)` (PseudoDojo publishes these
three levels; families with a single value give it for all three). They
come from PseudoPotentialData in the DFTK environment, read by a short
helper process (which downloads the family the first time), and are kept
for the session. Elements the family lacks, or a family without hints
(e.g. a `{El}` file pattern), are left out.
"""
function dftk_cutoff_hints(elements; family::AbstractString=DFTK_PSEUDO_FAMILY,
                           env::AbstractString=default_dftk_env())
    els = unique(Symbol.(collect(elements)))
    occursin("{El}", family) && return Dict{Symbol,eltype(values(CUTOFF_HINTS))}()
    missing_els = filter(e -> !haskey(CUTOFF_HINTS, (family, e)), els)
    if !isempty(missing_els) && dftk_env_ready(env)
        code = """
        using PseudoPotentialData, TOML
        f = PseudoFamily(ARGS[1]); out = Dict{String,Any}()
        for el in ARGS[2:end]
            m = try pseudometa(f, Symbol(el)) catch; continue end
            e0 = get(m, "Ecut", 0)
            out[el] = Dict(k => get(get(m, "cutoffs_" * k, Dict()), "Ecut", e0) for k in ("low", "normal", "high"))
        end
        TOML.print(out)
        """
        text = try
            read(pipeline(julia_in(env, "-e", code, family, string.(missing_els)...); stderr=devnull), String)
        catch
            ""
        end
        merge!(CUTOFF_HINTS, parse_cutoff_hints(text, family))
    end
    return Dict(e => CUTOFF_HINTS[(family, e)] for e in els if haskey(CUTOFF_HINTS, (family, e)))
end

"""Hints from the helper's TOML output (element tables with low/normal/high Ecut, Ha)."""
function parse_cutoff_hints(text::AbstractString, family::AbstractString)
    d = try
        TOML.parse(text)
    catch
        Dict{String,Any}()
    end
    out = Dict{Tuple{String,Symbol},NamedTuple{(:low, :normal, :high),NTuple{3,Float64}}}()
    for (el, v) in d
        v isa Dict && all(k -> get(v, k, 0) isa Real && get(v, k, 0) > 0, ("low", "normal", "high")) || continue
        out[(family, Symbol(el))] = (low=Float64(v["low"]), normal=Float64(v["normal"]), high=Float64(v["high"]))
    end
    return out
end

"""
    recommended_ecut(c; precision=:normal, family, env) -> cutoff in Ry, or nothing

The cutoff for `c`: the largest recommendation ([`dftk_cutoff_hints`](@ref))
among its elements at `precision` (`:low`, `:normal` or `:high`), in Ry
(twice the Ha value, as the tools and `ecutwfc` use Ry). `nothing` when
there are no hints for some element.
"""
function recommended_ecut(c::Crystal; precision::Symbol=:normal, family::AbstractString=DFTK_PSEUDO_FAMILY,
                          env::AbstractString=default_dftk_env())
    precision in (:low, :normal, :high) || throw(ArgumentError("precision must be :low, :normal or :high"))
    els = species_list(c)
    hints = dftk_cutoff_hints(els; family, env)
    all(e -> haskey(hints, e), els) || return nothing
    return 2 * maximum(getproperty(hints[e], precision) for e in els)
end

"""Stop a running job."""
function stop_job!(job::DFTKJob)
    process_running(job.process) || return job
    job.stopped = true
    kill(job.process)
    return job
end

"""Wait for the job to end; returns its [`job_status`](@ref)."""
function Base.wait(job::DFTKJob)
    wait(job.process)
    return job_status(job)
end

"""
The last `n` non-empty lines of the job's log (progress, or the error), as
plain text: terminal colour codes and other control characters removed.
"""
function job_log(job::DFTKJob; n::Integer=10)
    isfile(job.log) || return String[]
    lines = filter(!isempty ∘ strip, plain_text.(readlines(job.log)))
    return lines[max(1, end - n + 1):end]
end

"""`s` without ANSI escape sequences (colours) and other control characters."""
plain_text(s::AbstractString) = filter(!iscntrl, replace(replace(s, r"\e\[[0-9;?]*[ -/]*[@-~]" => ""), '\t' => ' '))

"""
    DFTKResults

What a finished DFTK run reports (from `<name>.results.toml`). Energies in
eV (`energy`, `terms`, `fermi`, `homo`, `lumo`; `nothing` where not known),
`forces` in eV/Å per atom, `structure` the final structure (the relaxed one
after a relaxation), `calculation` `"scf"` or `"relax"`, `converged`,
`iterations` (SCF iterations or relaxation steps) and, for an SCF,
`eigenvalues` (eV, per k-point and spin), `occupations` and `kweights`.
With the extra non-self-consistent step (`run_dftk(…; bands=true)`):
`dense`, the levels on a denser k-grid ([`DenseLevels`](@ref)), and `path`,
the band structure along the high-symmetry path ([`BandPath`](@ref));
`nothing` otherwise. After a relaxation, `history` holds every step
([`RelaxHistory`](@ref)); with `bands=true` the levels (and `dense`, `path`)
are those of an SCF of the relaxed structure.
"""
struct DFTKResults
    calculation::String
    converged::Bool
    iterations::Int
    energy::Float64
    terms::Vector{Pair{String,Float64}}
    fermi::Union{Nothing,Float64}
    homo::Union{Nothing,Float64}
    lumo::Union{Nothing,Float64}
    forces::Vector{Vec3}
    structure::Crystal
    eigenvalues::Vector{Vector{Float64}}
    occupations::Vector{Vector{Float64}}
    kweights::Vector{Float64}
    n_spin::Int
    version::String
    dense::Any                # Union{Nothing,DenseLevels}
    path::Any                 # Union{Nothing,BandPath}
    history::Any              # Union{Nothing,RelaxHistory}
end

"""
    RelaxHistory

The steps of a relaxation: `energies` (eV), the largest force on an atom
that may move (`max_forces`, eV/Å; clamped atoms keep their forces and are
left out), the structure of each step (`structures`) and the force on every
atom (`forces`, eV/Å), and which atoms were `clamped`.
"""
struct RelaxHistory
    energies::Vector{Float64}
    max_forces::Vector{Float64}
    structures::Vector{Crystal}
    forces::Vector{Vector{Vec3}}
    clamped::Vector{Bool}
end

Base.length(h::RelaxHistory) = length(h.energies)

"""The relaxation steps from `<name>.relax.toml` (species and periodicity from `like`), or `nothing`."""
function read_relax_history(path::AbstractString, like::Crystal)
    isfile(path) || return nothing
    d = TOML.parsefile(path)
    steps = get(d, "steps", Any[])
    isempty(steps) && return nothing
    n = natoms(like)
    clamped = Bool.(get(d, "clamped", fill(false, n)))
    structures = map(enumerate(steps)) do (k, s)
        L = Lattice((Float64.(v) for v in s["cell_A"])...)
        Crystal(L, like.species, [Vec3(Float64.(f)) for f in s["positions_frac"]]; pbc=like.pbc,
                title=like.title * " (DFTK step $(s["step"]))")
    end
    forces = [[Vec3(HA_PER_BOHR_EV_PER_A .* Float64.(f)) for f in s["forces_Ha_per_bohr"]] for s in steps]
    free = findall(!, clamped)
    maxf = [isempty(free) ? 0.0 : maximum(norm, f[free]) for f in forces]
    return RelaxHistory([HARTREE_EV * s["energy_Ha"] for s in steps], maxf, structures, forces, clamped)
end

"""
    DenseLevels

Levels (eV) of a non-self-consistent run on a denser Γ-centred k-grid
`kgrid`: `eigenvalues` per irreducible k-point (spin up first, then spin
down), with `kweights` summing to 1 per spin. [`dos`](@ref) uses them when
present.
"""
struct DenseLevels
    kgrid::NTuple{3,Int}
    eigenvalues::Vector{Vector{Float64}}
    kweights::Vector{Float64}
end

"""
    BandPath

Band structure along the high-symmetry path: `distances` along the path
(1/Å, 2π included), `branch_starts` (index where each continuous piece
starts; the path jumps between pieces), tick positions `ticks` with their
`labels` (e.g. "Γ", "X", "U | K" at a jump) and `eigenvalues[σ][k, band]`
(eV) per spin.
"""
struct BandPath
    distances::Vector{Float64}
    branch_starts::Vector{Int}
    ticks::Vector{Float64}
    labels::Vector{String}
    eigenvalues::Vector{Matrix{Float64}}
end

const HARTREE_EV = 27.211386245981          # eV (CODATA 2022)
const HA_PER_BOHR_EV_PER_A = HARTREE_EV / BOHR

"""
    read_dftk_results(path; title="") -> DFTKResults

Read the results file written by a DFTK script from [`write_dftk`](@ref)
with `results=…` (as [`run_dftk`](@ref) does), and the band file next to it
(`<name>.bands.toml`) when there is one. `homo` and `lumo` are the highest
SCF level below the Fermi level and the lowest above it.
"""
function read_dftk_results(path::AbstractString; title::AbstractString="",
                           bands::AbstractString=replace(path, r"\.results\.toml$" => ".bands.toml"),
                           history::AbstractString=replace(path, r"\.results\.toml$" => ".relax.toml"))
    d = TOML.parsefile(path)
    L = Lattice((Float64.(v) for v in d["cell_A"])...)
    pos = [Vec3(Float64.(p)) for p in d["positions_A"]]
    pbc = Tuple(Bool.(get(d, "periodic", [true, true, true])))
    c = crystal_from_cartesian(L, Symbol.(d["species"]), pos; pbc, title)
    ev(x) = x === nothing ? nothing : HARTREE_EV * x
    eigs = [HARTREE_EV .* Float64.(e) for e in get(d, "eigenvalues_Ha", Any[])]
    fermi = ev(get(d, "fermi_level_Ha", nothing))
    all_e = reduce(vcat, eigs; init=Float64[])
    below = fermi === nothing ? Float64[] : filter(<=(fermi), all_e)
    above = fermi === nothing ? Float64[] : filter(>(fermi), all_e)
    terms = [k => HARTREE_EV * Float64(x) for (k, x) in get(d, "energy_terms_Ha", Dict())]
    order = ["Kinetic", "AtomicLocal", "AtomicNonlocal", "Ewald", "PspCorrection", "Hartree", "Xc", "Entropy"]
    sort!(filter!(p -> p.first != "total", terms); by=p -> something(findfirst(==(p.first), order), 100))
    return DFTKResults(d["calculation"], d["converged"], get(d, "scf_iterations", get(d, "steps", 0)),
                       HARTREE_EV * d["energy_Ha"], terms, fermi,
                       isempty(below) ? nothing : maximum(below), isempty(above) ? nothing : minimum(above),
                       [Vec3(HA_PER_BOHR_EV_PER_A .* Float64.(f)) for f in d["forces_Ha_per_bohr"]], c,
                       eigs, [Float64.(o) for o in get(d, "occupations", Any[])],
                       Float64.(get(d, "kpoint_weights", Float64[])), get(d, "n_spin", 1), get(d, "version", ""),
                       read_bands_file(bands)..., read_relax_history(history, c))
end

"""The dense levels and the band path from a `<name>.bands.toml` file, or `nothing`s."""
function read_bands_file(path::AbstractString)
    isfile(path) || return (nothing, nothing)
    b = TOML.parsefile(path)
    dense = haskey(b, "dos_eigenvalues_Ha") ?
            DenseLevels(Tuple(Int.(b["dos_kgrid"])), [HARTREE_EV .* Float64.(e) for e in b["dos_eigenvalues_Ha"]],
                        Float64.(b["dos_kpoint_weights"])) : nothing
    kpath = haskey(b, "path_eigenvalues_Ha") ?
            BandPath(Float64.(b["path_distances_inv_A"]), Int.(b["path_branch_starts"]),
                     Float64.(b["path_tick_distances_inv_A"]), String.(b["path_tick_labels"]),
                     [HARTREE_EV .* permutedims(reduce(hcat, [Float64.(e) for e in σ])) for σ in b["path_eigenvalues_Ha"]]) :
            nothing
    return (dense, kpath)
end

"""The results of a finished job, or `nothing` while it runs or when it failed."""
function job_results(job::DFTKJob)
    job_status(job) === :done && isfile(job.results) || return nothing
    title = job.crystal === nothing ? "" : (isempty(job.crystal.title) ? formula(job.crystal) : job.crystal.title)
    return read_dftk_results(job.results; title)
end

"""
    dos(r::DFTKResults; width=0.1, step=width/10, range=nothing) -> (; energies, dos, fermi)

Density of states from the levels of an SCF, each broadened into a Gaussian
of standard deviation `width` (eV) and weighted with its k-point weight:
`dos[i, σ]` is in states per eV per cell for spin σ (one column without spin
polarization, where each level holds two electrons; spin up and down
otherwise), at `energies[i]` (eV, absolute; `fermi` is the Fermi level). By
default the energies cover all computed levels with 4 widths to spare.

The levels are those of the denser k-grid ([`DenseLevels`](@ref)) when the
run computed them (`levels=:dense`, the default `:auto`), which also include
more empty bands; otherwise (`:scf`) those of the SCF, whose few empty bands
make the DOS stop a little above the Fermi level.
"""
function dos(r::DFTKResults; width::Real=0.1, step::Real=width / 10, range=nothing, levels::Symbol=:auto)
    eigs, kw = dos_levels(r, levels)
    isempty(eigs) && throw(ArgumentError("no levels in these results (a relaxation keeps none)"))
    width > 0 && step > 0 || throw(ArgumentError("width and step must be positive"))
    nspin = r.n_spin
    nk = length(eigs) ÷ nspin
    w = length(kw) == length(eigs) ? kw : fill(1.0 / nk, length(eigs))
    lo, hi = range === nothing ? (minimum(minimum, eigs) - 4width, maximum(maximum, eigs) + 4width) :
             Float64.(range)
    E = collect(lo:step:hi)
    D = zeros(length(E), nspin)
    filled = 2 / nspin                                     # electrons per level and k-point
    g = 1 / (width * sqrt(2π))
    for σ in 1:nspin, ik in (σ - 1) * nk + 1:σ * nk, ε in eigs[ik]
        @. D[:, σ] += filled * w[ik] * g * exp(-((E - ε) / width)^2 / 2)
    end
    return (; energies=E, dos=D, fermi=r.fermi)
end

"""The levels and k-weights the DOS uses: `:auto` (dense when there), `:dense` or `:scf`."""
function dos_levels(r::DFTKResults, which::Symbol=:auto)
    which in (:auto, :dense, :scf) || throw(ArgumentError("levels must be :auto, :dense or :scf"))
    which === :dense && r.dense === nothing && throw(ArgumentError("no dense-grid levels in these results"))
    use_dense = which === :dense || (which === :auto && r.dense !== nothing)
    return use_dense ? (r.dense.eigenvalues, r.dense.kweights) : (r.eigenvalues, r.kweights)
end

"""
    write_dos(path, r; width=0.1) -> path

The DOS ([`dos`](@ref)) as CSV: energy (eV), energy relative to the Fermi
level, and the DOS (states/eV per cell) per spin.
"""
function write_dos(path::AbstractString, r::DFTKResults; kw...)
    d = dos(r; kw...)
    ef = something(d.fermi, 0.0)
    open(path, "w") do io
        println(io, "energy_eV,energy_minus_fermi_eV,", size(d.dos, 2) == 1 ? "dos_per_eV" : "dos_up_per_eV,dos_down_per_eV")
        for (i, e) in enumerate(d.energies)
            println(io, @sprintf("%.5f,%.5f,", e, e - ef), join((@sprintf("%.6f", x) for x in d.dos[i, :]), ","))
        end
    end
    return path
end

"""
    is_metallic(r; tol=0.003) -> Bool

Whether the levels make a metal: some band (the n-th level of each k-point,
per spin) lies above the Fermi level at one k-point and below it at another
(more than `tol` eV away; as DFTK's `is_metal`, over the SCF and the dense
k-grid), or an SCF level within `tol` of the Fermi level is partly filled
(e.g. a molecule's half-filled level). Partly filled levels alone are not
enough: with smearing, the Fermi level of an insulator can sit at the top of
the valence band and leave it slightly empty (GaAs at Γ).
"""
function is_metallic(r::DFTKResults; tol::Real=0.003)
    r.fermi === nothing && return false
    ef, full = r.fermi, 2 / r.n_spin
    for (eigs, _) in (dos_levels(r, :scf), r.dense === nothing ? (Vector{Float64}[], Float64[]) : dos_levels(r, :dense))
        isempty(eigs) && continue
        nk = length(eigs) ÷ r.n_spin
        for σ in 1:r.n_spin
            ks = (σ - 1) * nk + 1:σ * nk
            for b in 1:minimum(k -> length(eigs[k]), ks)
                any(k -> eigs[k][b] > ef + tol, ks) && any(k -> eigs[k][b] < ef - tol, ks) && return true
            end
        end
    end
    return any(zip(r.eigenvalues, r.occupations)) do (εs, os)
        any(abs(ε - ef) <= tol && 0.01full < o < 0.99full for (ε, o) in zip(εs, os))
    end
end

"""
    band_edges(r) -> (; vbm, cbm, gap, direct, vbm_at, cbm_at, source) or nothing

Valence-band maximum and conduction-band minimum (eV): the highest level
below the Fermi level and the lowest above it, over every k-point the run
computed (SCF grid, dense grid, band path). With a band path, `direct` says
whether both lie at the same point of the path and `vbm_at`/`cbm_at` where
("Γ", "between Γ and X", …); otherwise those are `nothing`. `source` says
which levels were available. `nothing` for a metal ([`is_metallic`](@ref))
or without levels above or below the Fermi level.
"""
function band_edges(r::DFTKResults)
    (r.fermi === nothing || is_metallic(r)) && return nothing
    ef = r.fermi
    sets = Vector{Float64}[collect(Iterators.flatten(r.eigenvalues))]
    r.dense === nothing || push!(sets, collect(Iterators.flatten(r.dense.eigenvalues)))
    p = r.path
    p === nothing || append!(sets, [vec(m) for m in p.eigenvalues])
    all_e = reduce(vcat, sets)
    below, above = filter(<=(ef), all_e), filter(>(ef), all_e)
    (isempty(below) || isempty(above)) && return nothing
    vbm, cbm = maximum(below), minimum(above)
    source = p !== nothing ? "band structure and k-grids" : r.dense !== nothing ? "dense k-grid" : "SCF k-grid"
    direct = vbm_at = cbm_at = nothing
    if p !== nothing
        kv = argmax(k -> maximum(m[k, b] for m in p.eigenvalues for b in axes(m, 2) if m[k, b] <= ef; init=-Inf),
                    eachindex(p.distances))
        kc = argmin(k -> minimum(m[k, b] for m in p.eigenvalues for b in axes(m, 2) if m[k, b] > ef; init=Inf),
                    eachindex(p.distances))
        vbm_at, cbm_at = path_point_name(p, kv), path_point_name(p, kc)
        direct = isapprox(p.distances[kv], p.distances[kc]; atol=1e-6) || vbm_at == cbm_at
    end
    return (; vbm, cbm, gap=cbm - vbm, direct, vbm_at, cbm_at, source)
end

"""The name of point `k` of the path: its label, or "between A and B"."""
function path_point_name(p::BandPath, k::Integer)
    x = p.distances[k]
    i = findfirst(t -> isapprox(t, x; atol=1e-6), p.ticks)
    i === nothing || return p.labels[i]
    j = something(findlast(<(x), p.ticks), 1)
    return "between " * split(p.labels[j], " | ")[end] * " and " * split(p.labels[min(j + 1, end)], " | ")[1]
end

"""
Band gap in eV ([`band_edges`](@ref): over all computed k-points), or
`nothing` when it is not known or the levels are partly filled
([`is_metallic`](@ref)). On a k-grid alone it can be larger than the true
gap; the band path finds the edges between grid points.
"""
band_gap(r::DFTKResults) = (e = band_edges(r); e === nothing ? nothing : e.gap)

"""Largest force on an atom (eV/Å), leaving out atoms clamped in a relaxation."""
function max_force(r::DFTKResults)
    free = r.history === nothing ? eachindex(r.forces) : findall(!, r.history.clamped)
    return isempty(r.forces) || isempty(free) ? 0.0 : maximum(norm, r.forces[free])
end

"""A few lines summarising the results."""
function summarize(r::DFTKResults)
    lines = [@sprintf("%s %s (%d %s)", r.calculation == "scf" ? "SCF" : "Relaxation",
                      r.converged ? "converged" : "NOT converged", r.iterations,
                      r.calculation == "scf" ? "iterations" : "steps"),
             @sprintf("Total energy %.6f eV (%.8f Ha)", r.energy, r.energy / HARTREE_EV),
             @sprintf("Largest force %.4f eV/Å", max_force(r))]
    h = r.history
    if h !== nothing && length(h) > 1
        push!(lines, @sprintf("Energy change %.4f eV over %d steps; largest force on a free atom %.3f → %.3f eV/Å",
                              last(h.energies) - first(h.energies), length(h) - 1, first(h.max_forces), last(h.max_forces)))
    end
    r.fermi === nothing || push!(lines, @sprintf("Fermi level %.4f eV", r.fermi))
    e = band_edges(r)
    if is_metallic(r)
        push!(lines, "Metallic: partly filled levels at the Fermi level")
    elseif e !== nothing && e.direct !== nothing
        push!(lines, @sprintf("Band gap %.4f eV, %s", e.gap, e.direct ? "direct at $(e.vbm_at)" :
                                                           "indirect: VBM $(e.vbm_at), CBM $(e.cbm_at)"))
    elseif e !== nothing
        push!(lines, @sprintf("Gap %.4f eV (on the %s)", e.gap, e.source))
    end
    return join(lines, "\n")
end

function Base.show(io::IO, ::MIME"text/plain", r::DFTKResults)
    println(io, "DFTK results for ", formula(r.structure), " (DFTK ", r.version, ")")
    print(io, summarize(r))
end

# --- convergence check -----------------------------------------------------------------

"""Cutoffs (Ry) for a convergence check: around the recommended one, or a fixed series without hints."""
function convergence_cutoffs(recommended::Union{Nothing,Real})
    recommended === nothing && return [20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0]
    return unique([2 * round(f * recommended / 2) for f in (0.5, 0.625, 0.75, 0.875, 1.0, 1.25, 1.5)])
end

"""k-spacings (1/Å) for a convergence check, coarse to fine, keeping one per distinct k-grid of `c`."""
function convergence_kspacings(c::Crystal; spacings=(0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15))
    seen, out = Set{NTuple{3,Int}}(), Float64[]
    for s in spacings
        g = kpoint_grid(c; kspacing=s)
        g in seen && continue
        push!(seen, g); push!(out, s)
    end
    return out
end

"""
    write_dftk_convergence(io, c; ecuts, kspacings, ecutwfc, kspacing, pseudo, results)

A DFTK script that runs one SCF per cutoff in `ecuts` (Ry) at the k-grid of
`kspacing`, then one per k-spacing in `kspacings` at the cutoff `ecutwfc`
(Ry), and writes energy, Fermi level, band edges and SCF convergence of each
to `results` (TOML) after every step, for [`read_convergence`](@ref).
Constraints are ignored (no relaxation).
"""
function write_dftk_convergence(io::IO, c::Crystal; ecuts, kspacings, ecutwfc::Real=40, kspacing::Real=0.25,
                                pseudo::AbstractString=DFTK_PSEUDO_FAMILY, results::AbstractString="convergence.toml")
    title = isempty(c.title) ? formula(c) : c.title
    println(io, "# ", title, " — DFTK.jl convergence check written by AtomBond")
    println(io, "using DFTK, PseudoPotentialData, AtomsBase, Unitful, UnitfulAtomic")
    println(io, "import TOML")
    write_dftk_system(io, c, pseudo)
    println(io, "model = model_DFT(system; model_kwargs...)")
    println(io, "\n# (series, cutoff in Ry, k-spacing in 1/Å, k-grid)")
    println(io, "steps = [")
    for e in ecuts
        @printf(io, "    (\"cutoff\", %.2f, %.4f, [%d, %d, %d]),\n", e, kspacing, kpoint_grid(c; kspacing)...)
    end
    for s in kspacings
        @printf(io, "    (\"kpoints\", %.2f, %.4f, [%d, %d, %d]),\n", ecutwfc, s, kpoint_grid(c; kspacing=s)...)
    end
    println(io, "]")
    print(io, """
    rows = Dict{String,Any}[]
    for (series, ecut_Ry, kspacing, kgrid) in steps
        basis = PlaneWaveBasis(model; Ecut=ecut_Ry / 2, kgrid)
        scfres = self_consistent_field(basis; ρ=guess_density(basis, system), tol=1e-8)
        eigs = reduce(vcat, scfres.eigenvalues)
        below, above = filter(<=(scfres.εF), eigs), filter(>(scfres.εF), eigs)
        row = Dict{String,Any}("series" => series, "ecut_Ry" => ecut_Ry, "kspacing" => kspacing, "kgrid" => kgrid,
                               "n_kpoints" => length(basis.kpoints), "energy_Ha" => scfres.energies.total,
                               "converged" => scfres.converged, "iterations" => scfres.n_iter,
                               "fermi_level_Ha" => scfres.εF)
        isempty(below) || (row["homo_Ha"] = maximum(below))
        isempty(above) || (row["lumo_Ha"] = minimum(above))
        push!(rows, row)
        out = Dict{String,Any}("program" => "DFTK", "version" => string(pkgversion(DFTK)), "natoms" => length(system),
                               "planned" => length(steps), "steps" => rows)
        open(io -> TOML.print(io, out), joinpath(@__DIR__, $(repr(String(results)))), "w")
        println("AtomBond convergence step ", length(rows), "/", length(steps), ": ", series, " ",
                ecut_Ry, " Ry, k-grid ", join(kgrid, "×"), ", E = ", scfres.energies.total, " Ha")
    end
    """)
end

"""
    run_dftk_convergence(c; ecuts=auto, kspacings=auto, ecutwfc=auto, kspacing=0.25, dir, name, env, pseudo) -> DFTKJob

Start a convergence check of `c` in DFTK (script from `write_dftk_convergence`):
SCF total energies for a series of cutoffs (by default around the cutoff the
pseudopotentials recommend, [`recommended_ecut`](@ref), from half to 1.5
times it) at the k-grid of `kspacing`, and for a series of k-spacings (one
per distinct k-grid; none for a molecule) at the cutoff `ecutwfc` (default:
the recommended one, or 40 Ry). The results go to `name.convergence.toml`;
read them with [`job_convergence`](@ref) or [`read_convergence`](@ref).
"""
function run_dftk_convergence(c::Crystal; name::AbstractString=formula(c),
                              dir::AbstractString=joinpath(pwd(), "dftk_" * name * "_conv_" *
                                                           Dates.format(Dates.now(), "yyyymmdd-HHMMSS")),
                              env::AbstractString=default_dftk_env(), pseudo::AbstractString=DFTK_PSEUDO_FAMILY,
                              ecuts=nothing, kspacings=nothing, ecutwfc=nothing, kspacing::Real=0.25)
    dftk_env_ready(env) || throw(ArgumentError("no DFTK environment at $env; create it with install_dftk()"))
    occursin(r"^[\w.+-]+$", name) || throw(ArgumentError("use letters, digits, . + - _ in the name"))
    rec = ecuts === nothing || ecutwfc === nothing ? recommended_ecut(c; family=pseudo, env) : nothing
    ecuts = ecuts === nothing ? convergence_cutoffs(rec) : Float64.(collect(ecuts))
    ecutwfc = ecutwfc === nothing ? something(rec, 40.0) : Float64(ecutwfc)
    kspacings = kspacings === nothing ? convergence_kspacings(c) : Float64.(collect(kspacings))
    length(unique(kpoint_grid(c; kspacing=s) for s in kspacings)) > 1 || (kspacings = Float64[])
    mkpath(dir)
    script, log, results = (joinpath(dir, name * ext) for ext in (".conv.jl", ".conv.log", ".convergence.toml"))
    isfile(results) && rm(results)
    open(io -> write_dftk_convergence(io, c; ecuts, kspacings, ecutwfc, kspacing, pseudo, results=basename(results)),
         script, "w")
    proc = open(log, "w") do io
        run(pipeline(julia_in(env, script; dir); stdout=io, stderr=io); wait=false)
    end
    return DFTKJob(c, dir, script, log, results, proc, time(), false)
end

"""One SCF of a convergence check: cutoff (Ry), k-spacing (1/Å), k-grid, energy (eV), gap (eV) or nothing."""
struct ConvergenceStep
    series::String
    ecut::Float64
    kspacing::Float64
    kgrid::NTuple{3,Int}
    nkpoints::Int
    energy::Float64
    converged::Bool
    gap::Union{Nothing,Float64}
end

"""
    ConvergenceResults

A convergence check read by [`read_convergence`](@ref): `natoms`, the SCF
`steps` ([`ConvergenceStep`](@ref)) of the `"cutoff"` and `"kpoints"`
series, `planned` (all steps of the run; fewer are present while it runs)
and the DFTK `version`.
"""
struct ConvergenceResults
    natoms::Int
    steps::Vector{ConvergenceStep}
    planned::Int
    version::String
end

"""Read a `name.convergence.toml` (also while the check is still running)."""
function read_convergence(path::AbstractString)
    d = TOML.parsefile(path)
    steps = map(get(d, "steps", Any[])) do s
        gap = haskey(s, "homo_Ha") && haskey(s, "lumo_Ha") ? HARTREE_EV * (s["lumo_Ha"] - s["homo_Ha"]) : nothing
        ConvergenceStep(s["series"], s["ecut_Ry"], s["kspacing"], Tuple(Int.(s["kgrid"])), get(s, "n_kpoints", 0),
                        HARTREE_EV * s["energy_Ha"], s["converged"], gap)
    end
    return ConvergenceResults(d["natoms"], steps, get(d, "planned", length(steps)), get(d, "version", ""))
end

"""The convergence results of a job so far (`nothing` before the first step)."""
job_convergence(job::DFTKJob) = isfile(job.results) && endswith(job.results, ".convergence.toml") ?
                                read_convergence(job.results) : nothing

"""
    convergence_series(cr, series) -> (; parameter, energy_error, gap, steps)

The steps of one series (`"cutoff"`: by cutoff, ascending; `"kpoints"`: by
k-spacing, coarse to fine) with each energy's difference from the most
accurate one (highest cutoff, densest grid) in meV per atom.
"""
function convergence_series(cr::ConvergenceResults, series::AbstractString)
    st = filter(s -> s.series == series, cr.steps)
    sort!(st; by=s -> series == "cutoff" ? s.ecut : -s.kspacing)
    isempty(st) && return (; parameter=Float64[], energy_error=Float64[], gap=Union{Nothing,Float64}[], steps=st)
    ref = last(st).energy
    return (; parameter=[series == "cutoff" ? s.ecut : s.kspacing for s in st],
            energy_error=[1000 * (s.energy - ref) / cr.natoms for s in st], gap=[s.gap for s in st], steps=st)
end

"""
    converged_settings(cr; tol=1.0) -> (; ecut, kspacing)

The smallest cutoff (Ry) and the largest k-spacing (1/Å) whose total energy,
and that of every more accurate step of the series, is within `tol` meV per
atom of the most accurate one; `nothing` for a series that is missing or
not converged within its range (only its most accurate step qualifies).
"""
function converged_settings(cr::ConvergenceResults; tol::Real=1.0)
    function first_ok(series)
        s = convergence_series(cr, series)
        length(s.steps) < 2 && return nothing
        for i in 1:length(s.steps) - 1
            all(abs(e) <= tol for e in s.energy_error[i:end]) && return s.parameter[i]
        end
        return nothing
    end
    return (; ecut=first_ok("cutoff"), kspacing=first_ok("kpoints"))
end

"""A few lines summarising a convergence check."""
function summarize(cr::ConvergenceResults; tol::Real=1.0)
    rec = converged_settings(cr; tol)
    lines = ["Convergence check: $(length(cr.steps)) of $(cr.planned) SCF runs" *
             (any(s -> !s.converged, cr.steps) ? " (some SCF NOT converged)" : "")]
    for (series, what, unit) in (("cutoff", "Cutoff", "Ry"), ("kpoints", "k-spacing", "1/Å"))
        s = convergence_series(cr, series)
        isempty(s.steps) && continue
        val = series == "cutoff" ? rec.ecut : rec.kspacing
        push!(lines, val === nothing ? @sprintf("%s: not within %g meV/atom below %g %s", what, tol, last(s.parameter), unit) :
                     @sprintf("%s: %g %s is within %g meV/atom", what, val, unit, tol) *
                     (series == "kpoints" ? " (k-grid " * join(s.steps[findfirst(==(val), s.parameter)].kgrid, "×") * ")" : ""))
    end
    return join(lines, "\n")
end

function Base.show(io::IO, ::MIME"text/plain", cr::ConvergenceResults)
    println(io, "DFTK convergence check (DFTK ", cr.version, ")")
    print(io, summarize(cr))
end
