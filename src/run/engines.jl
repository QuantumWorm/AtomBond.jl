# DFT engines: the codes AtomBond can run (DFTK by default, Quantum ESPRESSO,
# CASTEP, and codes users add themselves). AtomBond holds only the adapters
# (input writer, command, output reader); the programs and their
# pseudopotentials stay on the user's machine, found through a settings file,
# an environment variable or the PATH (see `find_engine`).

"""
    DFTEngine

A code AtomBond can run. An engine is a (usually empty) struct with methods
for [`engine_id`](@ref), [`engine_name`](@ref),
[`write_engine_inputs`](@ref) and [`read_engine_results`](@ref), and
optionally [`engine_executables`](@ref), [`engine_command`](@ref),
[`engine_tasks`](@ref), [`engine_finished`](@ref),
[`engine_history`](@ref), [`engine_pseudo_default`](@ref) and
[`engine_install_hint`](@ref); [`register_engine!`](@ref) adds it to the Run
tool. See *Adding a DFT engine* in the manual for a template.
"""
abstract type DFTEngine end

"""The engines in the order the Run tool lists them (DFTK first)."""
const DFT_ENGINES = DFTEngine[]

"""
    register_engine!(engine) -> engine

Add `engine` to the engines AtomBond knows (replacing one with the same
[`engine_id`](@ref)), so that [`run_dft`](@ref) and the Run tool can use it.
"""
function register_engine!(e::DFTEngine)
    k = findfirst(x -> engine_id(x) === engine_id(e), DFT_ENGINES)
    k === nothing ? push!(DFT_ENGINES, e) : (DFT_ENGINES[k] = e)
    return e
end

"""`dft_engines()` lists the registered engines; `dft_engine(:qe)` gets one by its id."""
dft_engines() = copy(DFT_ENGINES)
function dft_engine(id::Symbol)
    k = findfirst(x -> engine_id(x) === id, DFT_ENGINES)
    k === nothing && throw(ArgumentError("no DFT engine :$id (known: " * join((":" * string(engine_id(e)) for e in DFT_ENGINES), ", ") * ")"))
    return DFT_ENGINES[k]
end
dft_engine(e::DFTEngine) = e

"""`engine_id(e)`: a short Symbol (`:qe`); used in the settings file and in `ATOMBOND_<ID>_COMMAND`."""
function engine_id end
"""`engine_name(e)`: the name shown to people ("Quantum ESPRESSO")."""
function engine_name end

"""
    write_engine_inputs(e, dir, name, c; relax, kspacing, ecutwfc, pseudo) -> (input, output)

Write the input files of a calculation of `c` into `dir` (they exist). `relax`
is `:none` (an SCF), `:atoms` or `:cell`; `ecutwfc` is in Ry, `kspacing` in
1/Å (2π included, see [`kpoint_grid`](@ref)); `pseudo` is the text of the
tool's *Pseudopotentials* field. Returns the main input file and the file
that exists once the calculation has finished (its results).
"""
function write_engine_inputs end

"""`read_engine_results(e, job) -> DFTResults`: the results of a finished job."""
function read_engine_results end

"""Program names to look for on the PATH, most preferred first."""
engine_executables(::DFTEngine) = String[]

"""
    engine_command(e, command, dir, name, input) -> Cmd

The command that runs the calculation (in `dir`, output to the job's log):
by default `command` followed by the input file name.
"""
engine_command(::DFTEngine, command::Vector{String}, dir, name, input) = `$command $(basename(input))`

"""Tasks of the Run tool the engine can do: some of `:scf`, `:bands`, `:relax_atoms`, `:relax_cell`, `:convergence`."""
engine_tasks(::DFTEngine) = [:scf, :relax_atoms, :relax_cell]

"""Whether a job whose process ended successfully left its results (by default: the output file exists)."""
engine_finished(::DFTEngine, job) = isfile(job.results)

"""The steps of a relaxation so far ([`RelaxHistory`](@ref)), or `nothing` (the default)."""
engine_history(::DFTEngine, job) = nothing

"""The default text of the tool's *Pseudopotentials* field for this engine."""
engine_pseudo_default(::DFTEngine) = ""

"""One line on how to get the program (shown when it is not found)."""
engine_install_hint(::DFTEngine) = ""

Base.show(io::IO, e::DFTEngine) = print(io, engine_name(e), " engine (:", engine_id(e), ")")

# --- where the programs are -------------------------------------------------------

"""
    engine_settings_file() -> String

The settings file of the engines: `\$ATOMBOND_ENGINES` when set, otherwise
`atombond/engines.toml` in `\$XDG_CONFIG_HOME` (default `~/.config`). One
table per engine id:

```toml
[castep]
command = "mpirun -np 4 castep.mpi"

[qe]
command = "/opt/qe/bin/pw.x"
pseudo_dir = "~/pseudo/sssp"
```
"""
function engine_settings_file()
    haskey(ENV, "ATOMBOND_ENGINES") && return ENV["ATOMBOND_ENGINES"]
    base = get(ENV, "XDG_CONFIG_HOME", joinpath(homedir(), ".config"))
    return joinpath(base, "atombond", "engines.toml")
end

"""All engine settings (an empty Dict without a settings file; a broken file is reported and ignored)."""
function engine_settings()
    f = engine_settings_file()
    isfile(f) || return Dict{String,Any}()
    return try
        TOML.parsefile(f)
    catch err
        @warn "Ignoring the engine settings in $f: $(sprint(showerror, err))"
        Dict{String,Any}()
    end
end

"""`engine_setting(e, key, default)`: one setting of an engine from the settings file."""
engine_setting(e, key::AbstractString, default=nothing) =
    get(get(engine_settings(), string(engine_id(dft_engine(e))), Dict{String,Any}()), key, default)

"""
    set_engine_setting!(e, key, value)

Store a setting of an engine (e.g. `"command"`, `"pseudo_dir"`) in the
settings file, keeping the others; `nothing` removes it.
"""
function set_engine_setting!(e, key::AbstractString, value)
    all = engine_settings()
    t = get!(all, string(engine_id(dft_engine(e))), Dict{String,Any}())
    value === nothing ? delete!(t, key) : (t[key] = value)
    f = engine_settings_file()
    mkpath(dirname(f))
    open(io -> TOML.print(io, all), f, "w")
    return value
end

"""
    find_engine(e) -> (; command, source, found) or nothing

The command that runs engine `e`, from the first of: the `command` in the
settings file ([`engine_settings_file`](@ref)), the environment variable
`ATOMBOND_<ID>_COMMAND` (e.g. `ATOMBOND_CASTEP_COMMAND="mpirun -np 4
castep.mpi"`), or the PATH ([`engine_executables`](@ref), e.g. `castep` in
`~/bin`). `command` is split into words like a shell would; `found` tells
whether its program exists. `nothing` when none of these gives a command.
"""
function find_engine(e)
    e = dft_engine(e)
    candidates = Tuple{String,String}[]
    s = engine_setting(e, "command")
    s isa AbstractString && !isempty(strip(s)) && push!(candidates, (s, "settings file"))
    var = "ATOMBOND_" * uppercase(string(engine_id(e))) * "_COMMAND"
    haskey(ENV, var) && !isempty(strip(ENV[var])) && push!(candidates, (ENV[var], var))
    isempty(candidates) || return command_info(candidates[1]...)
    for exe in engine_executables(e)
        p = Sys.which(exe)
        p === nothing || return (; command=[p], source="PATH", found=true)
    end
    return nothing
end

function command_info(text::AbstractString, source::AbstractString)
    words = [expanduser(w) for w in Base.shell_split(text)]
    # the program is the first word that is not an MPI launcher or its option
    prog = findfirst(w -> !startswith(w, "-") && !(w in ("mpirun", "mpiexec", "srun")) && tryparse(Int, w) === nothing, words)
    found = !isempty(words) && all(w -> Sys.which(w) !== nothing || isfile(w),
                                   (words[1], prog === nothing ? words[1] : words[prog]))
    return (; command=words, source, found)
end

"""One line saying where engine `e` was found, or that it was not."""
engine_status(id::Symbol; kw...) = engine_status(dft_engine(id); kw...)
function engine_status(e::DFTEngine)
    f = find_engine(e)
    if f === nothing
        hint = engine_install_hint(e)
        return "$(engine_name(e)) not found: set its command" * (isempty(hint) ? "" : "\n" * hint)
    end
    shown = replace(join(f.command, " "), homedir() => "~")
    return (f.found ? "" : "⚠ not found: ") * shown * "\n(from the " * (f.source == "PATH" ? "PATH" : f.source) * ")"
end

# --- running ------------------------------------------------------------------------

"""Resolve `relax=:auto` for `c`: `:none` without constraints, `:atoms` for a fixed cell, `:cell` otherwise."""
relax_mode(c::Crystal, relax::Symbol) =
    relax !== :auto ? relax : !has_constraints(c) ? :none : c.cell_constraints == (0, 0, 0, 0, 0, 0) ? :atoms : :cell

"""
    run_dft(c; engine=:dftk, dir, name, kspacing=0.25, ecutwfc=40, pseudo, relax=:auto, bands=false, command) -> DFTJob

Start a calculation of `c` with a DFT engine in the background ([`dft_engines`](@ref)):
writes the inputs into `dir` (default: a new folder `<engine>_<name>_<time>`),
runs the program there with its output in `name.log`, and returns the job
([`job_status`](@ref), [`job_results`](@ref), [`stop_job!`](@ref)). `relax`
as in [`run_dftk`](@ref): `:none`, `:atoms`, `:cell`, or `:auto` (relax when
`c` has constraints). `command` overrides the one [`find_engine`](@ref)
finds. With `engine=:dftk` this is [`run_dftk`](@ref) (which also does
`bands`).
"""
function run_dft(c::Crystal; engine=:dftk, name::AbstractString=formula(c),
                 dir::Union{Nothing,AbstractString}=nothing, kspacing::Real=0.25, ecutwfc::Real=40,
                 pseudo::Union{Nothing,AbstractString}=nothing, relax::Symbol=:auto, bands::Bool=false,
                 command=nothing, kw...)
    e = dft_engine(engine)
    relax in (:auto, :none, :atoms, :cell) || throw(ArgumentError("relax is :auto, :none, :atoms or :cell, got :$relax"))
    stamp = Dates.format(Dates.now(), "yyyymmdd-HHMMSS")
    folder = something(dir, joinpath(pwd(), string(engine_id(e)) * "_" * name * "_" * stamp))
    if engine_id(e) === :dftk
        return run_dftk(c; name, dir=folder, kspacing, ecutwfc, pseudo=something(pseudo, DFTK_PSEUDO_FAMILY), relax,
                        bands, kw...)
    end
    bands && throw(ArgumentError("band structures are only computed with DFTK so far"))
    occursin(r"^[\w.+-]+$", name) || throw(ArgumentError("use letters, digits, . + - _ in the name"))
    cmd = if command === nothing
        f = find_engine(e)
        f === nothing && throw(ArgumentError(engine_status(e)))
        f.command
    else
        command isa AbstractString ? command_info(command, "").command : collect(String, command)
    end
    mode = relax_mode(c, relax)
    mkpath(folder)
    input, output = write_engine_inputs(e, folder, name, c; relax=mode, kspacing, ecutwfc,
                                        pseudo=something(pseudo, engine_pseudo_default(e)))
    isfile(output) && rm(output)
    log = joinpath(folder, name * ".log")
    proc = open(log, "w") do io
        run(pipeline(Cmd(engine_command(e, cmd, folder, name, input); dir=folder); stdout=io, stderr=io); wait=false)
    end
    return DFTJob(c, folder, input, log, output, proc, time(), false, mode !== :none, engine_id(e))
end

# --- DFTK as an engine ------------------------------------------------------------

"""DFTK.jl, run in its own Julia environment ([`run_dftk`](@ref)); the default engine."""
struct DFTKEngine <: DFTEngine end
engine_id(::DFTKEngine) = :dftk
engine_name(::DFTKEngine) = "DFTK"
engine_tasks(::DFTKEngine) = [:bands, :scf, :relax_atoms, :relax_cell, :convergence]
engine_pseudo_default(::DFTKEngine) = DFTK_PSEUDO_FAMILY
engine_install_hint(::DFTKEngine) = "Run installs it (minutes, once)"
engine_finished(::DFTKEngine, job) = isfile(job.results)
function engine_status(e::DFTKEngine; env::AbstractString=default_dftk_env())
    return dftk_env_ready(env) ? "DFTK is installed there" : "No DFTK there yet: Run installs\nit first (minutes, once)"
end
function read_engine_results(::DFTKEngine, job)
    title = job.crystal === nothing ? "" : (isempty(job.crystal.title) ? formula(job.crystal) : job.crystal.title)
    return read_dftk_results(job.results; title)
end

register_engine!(DFTKEngine())
