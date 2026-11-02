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

"""A Julia command for a helper process: this Julia, no startup file, `env` as the project."""
function julia_in(env::AbstractString, args...; dir::AbstractString=pwd())
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$env $args`
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
"""
function run_dftk(c::Crystal; name::AbstractString=formula(c),
                  dir::AbstractString=joinpath(pwd(), "dftk_" * name * "_" * Dates.format(Dates.now(), "yyyymmdd-HHMMSS")),
                  env::AbstractString=default_dftk_env(), kspacing::Real=0.25, ecutwfc::Real=40,
                  pseudo::AbstractString=DFTK_PSEUDO_FAMILY)
    dftk_env_ready(env) || throw(ArgumentError("no DFTK environment at $env; create it with install_dftk()"))
    occursin(r"^[\w.+-]+$", name) || throw(ArgumentError("use letters, digits, . + - _ in the name"))
    mkpath(dir)
    script, log, results = (joinpath(dir, name * ext) for ext in (".jl", ".log", ".results.toml"))
    isfile(results) && rm(results)
    write_dftk(script, c; kspacing, ecutwfc, pseudo, results=basename(results))
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

"""The last `n` non-empty lines of the job's log (progress, or the error)."""
function job_log(job::DFTKJob; n::Integer=10)
    isfile(job.log) || return String[]
    lines = filter(!isempty ∘ strip, readlines(job.log))
    return lines[max(1, end - n + 1):end]
end

"""
    DFTKResults

What a finished DFTK run reports (from `<name>.results.toml`). Energies in
eV (`energy`, `terms`, `fermi`, `homo`, `lumo`; `nothing` where not known),
`forces` in eV/Å per atom, `structure` the final structure (the relaxed one
after a relaxation), `calculation` `"scf"` or `"relax"`, `converged`,
`iterations` (SCF iterations or relaxation steps) and, for an SCF,
`eigenvalues` (eV, per k-point and spin), `occupations` and `kweights`.
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
end

const HARTREE_EV = 27.211386245981          # eV (CODATA 2022)
const HA_PER_BOHR_EV_PER_A = HARTREE_EV / BOHR

"""
    read_dftk_results(path; title="") -> DFTKResults

Read the results file written by a DFTK script from [`write_dftk`](@ref)
with `results=…` (as [`run_dftk`](@ref) does). The band edges are the
highest eigenvalue below the Fermi level and the lowest above it.
"""
function read_dftk_results(path::AbstractString; title::AbstractString="")
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
                       Float64.(get(d, "kpoint_weights", Float64[])), get(d, "n_spin", 1), get(d, "version", ""))
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
default the energies cover all computed levels with 4 widths to spare; only
the bands DFTK computed are included, so the DOS stops a little above the
Fermi level.
"""
function dos(r::DFTKResults; width::Real=0.1, step::Real=width / 10, range=nothing)
    isempty(r.eigenvalues) && throw(ArgumentError("no levels in these results (a relaxation keeps none)"))
    width > 0 && step > 0 || throw(ArgumentError("width and step must be positive"))
    nspin = r.n_spin
    nk = length(r.eigenvalues) ÷ nspin
    w = length(r.kweights) == length(r.eigenvalues) ? r.kweights : fill(1.0 / nk, length(r.eigenvalues))
    lo, hi = range === nothing ? (minimum(minimum, r.eigenvalues) - 4width, maximum(maximum, r.eigenvalues) + 4width) :
             Float64.(range)
    E = collect(lo:step:hi)
    D = zeros(length(E), nspin)
    filled = 2 / nspin                                     # electrons per level and k-point
    g = 1 / (width * sqrt(2π))
    for σ in 1:nspin, ik in (σ - 1) * nk + 1:σ * nk, ε in r.eigenvalues[ik]
        @. D[:, σ] += filled * w[ik] * g * exp(-((E - ε) / width)^2 / 2)
    end
    return (; energies=E, dos=D, fermi=r.fermi)
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
Whether the SCF found partly filled levels (occupation between 1 % and 99 %
of a full level): a metal, or a system whose levels cross the Fermi level.
"""
function is_metallic(r::DFTKResults)
    full = 2 / r.n_spin
    return any(o -> 0.01full < o < 0.99full, Iterators.flatten(r.occupations))
end

"""
Band gap in eV (lowest level above the Fermi level minus highest below), or
`nothing` when it is not known or the levels are partly filled
([`is_metallic`](@ref)). Only the k-points of the grid count, so it can be
larger than the true gap.
"""
band_gap(r::DFTKResults) = r.homo === nothing || r.lumo === nothing || is_metallic(r) ? nothing : r.lumo - r.homo

"""Largest force on any atom (eV/Å)."""
max_force(r::DFTKResults) = isempty(r.forces) ? 0.0 : maximum(norm, r.forces)

"""A few lines summarising the results."""
function summarize(r::DFTKResults)
    lines = [@sprintf("%s %s (%d %s)", r.calculation == "scf" ? "SCF" : "Relaxation",
                      r.converged ? "converged" : "NOT converged", r.iterations,
                      r.calculation == "scf" ? "iterations" : "steps"),
             @sprintf("Total energy %.6f eV (%.8f Ha)", r.energy, r.energy / HARTREE_EV),
             @sprintf("Largest force %.4f eV/Å", max_force(r))]
    r.fermi === nothing || push!(lines, @sprintf("Fermi level %.4f eV", r.fermi))
    gap = band_gap(r)
    if is_metallic(r)
        push!(lines, "Metallic: partly filled levels at the Fermi level")
    elseif gap !== nothing
        push!(lines, @sprintf("Gap %.4f eV (on the k-point grid)", gap))
    end
    return join(lines, "\n")
end

function Base.show(io::IO, ::MIME"text/plain", r::DFTKResults)
    println(io, "DFTK results for ", formula(r.structure), " (DFTK ", r.version, ")")
    print(io, summarize(r))
end
