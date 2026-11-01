# Run DFTK tool: start a DFTK calculation of the current structure, follow its
# progress and show the results when it is done (src/run/dftk.jl does the work).

"""The calculation `write_dftk` will set up for `c`, in words."""
function dftk_plan(c::Crystal, kspacing)
    what = has_constraints(c) ? "Relaxation (the structure has constraints)" : "SCF (energies, forces, levels)"
    k = kspacing === nothing ? "" : ", k-grid " * join(kpoint_grid(c; kspacing), " × ")
    return what * "\n" * plural(natoms(c), "atom") * k
end

"""
    dftk_dialog(viewer) -> NamedTuple of the tool's controls

Run DFTK tool: runs the current structure in DFTK.jl in a separate Julia
process (see [`run_dftk`](@ref)) — an SCF, or a relaxation when the
structure has constraints. It shows the progress from the log while it runs
and a summary of the results when it is done; *Open result structure* opens
the final (relaxed) structure in a new tab. *Install DFTK* creates the DFTK
environment the first time ([`install_dftk`](@ref)).
"""
function dftk_dialog(v::Viewer)
    f = new_dialog(v, :dftk, "Run DFTK")
    f === nothing && return nothing
    heading!(f, "DFTK.jl calculation")
    plan = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], plan; halign=:left, fontsize=13, color=:gray25, justification=:left)

    heading!(f, "Environment")
    env = field!(f, "Julia environment", default_dftk_env(), s -> !isempty(strip(s)))
    envstate = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], envstate; halign=:left, fontsize=13, color=:gray25, justification=:left)
    install = button!(f, "Install DFTK")

    heading!(f, "Calculation")
    folder = field!(f, "Folder", pwd(), s -> !isempty(strip(s)))
    name = field!(f, "Name", formula(v.crystal), s -> occursin(r"^[\w.+-]+$", strip(s)))
    kspacing = field!(f, "k-spacing (1/Å, incl. 2π)", "0.25", s -> isnum(s) && parse(Float64, s) > 0)
    ecut = field!(f, "Cutoff (Ry)", "40", s -> isnum(s) && parse(Float64, s) > 0)
    pseudo = field!(f, "Pseudopotentials", DFTK_PSEUDO_FAMILY, s -> occursin("{El}", s) || occursin(r"^[\w.+-]+$", s))
    note!(f, "PBE, 1 mHa Gaussian smearing; cutoff/2 in Ha for DFTK.\nThe run goes into Folder/dftk_Name_<date-time>/.")
    r = nextrow!(f)
    run, stop = (with_click_feedback!(Makie.Button(f.g[r, k]; label, width=Makie.Relative(1), fontsize=PANEL_FONTSIZE,
                                                   padding=(6, 6, 5, 5))) for (k, label) in ((1, "Run"), (2, "Stop")))

    heading!(f, "Progress and results")
    progress = Makie.Observable("Not started")
    Makie.Label(f.g[nextrow!(f), 1:2], progress; halign=:left, fontsize=12, color=:gray25, justification=:left)
    summary = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], summary; halign=:left, fontsize=13, justification=:left)
    open_result = button!(f, "Open result structure")

    num(t) = (x = tryparse(Float64, strip(shown_text(t))); x === nothing || x <= 0 ? nothing : x)
    update_plan(_...) = (plan[] = dftk_plan(v.crystal, num(kspacing)))
    Makie.on(update_plan, kspacing.displayed_string)
    Makie.on(update_plan, v.header)                  # the structure changed
    update_plan()
    update_env(_...) = (envstate[] = dftk_env_ready(strip(shown_text(env))) ? "DFTK is installed there" :
                                     "No DFTK there yet: Install DFTK (several minutes, once)")
    Makie.on(update_env, env.displayed_string)
    update_env()

    state = get!(v.dialogs, :dftk_state, Dict{Symbol,Any}())    # the job outlives the tool
    controls = (; plan, env, envstate, install, folder, name, kspacing, ecut, pseudo, run, stop, progress, summary,
                open_result, state)
    show_job!(v, controls)

    Makie.on(install.clicks) do _
        e = form_value(v, env, "an environment folder")
        e === nothing && return
        busy(state) && (v.status[] = "Wait for the running job to end, or stop it"; return)
        dftk_env_ready(e) && (v.status[] = "DFTK is already installed in $e"; return)
        state[:job] = install_dftk(e)
        state[:results] = nothing
        v.status[] = "Installing DFTK into $e (log: install.log there)"
        follow_job!(v, controls)
    end
    Makie.on(run.clicks) do _
        vals = (form_value(v, env, "an environment folder"), form_value(v, folder, "a folder"),
                form_value(v, name, "a name"), form_value(v, kspacing, "a k-point spacing"),
                form_value(v, ecut, "a cutoff"), form_value(v, pseudo, "pseudopotentials"))
        any(isnothing, vals) && return
        busy(state) && (v.status[] = "A DFTK job is still running: stop it first"; return)
        e, dir, nm, ks, ec, ps = vals
        try
            stamp = Dates.format(Dates.now(), "yyyymmdd-HHMMSS")
            state[:job] = run_dftk(v.crystal; env=e, name=nm, dir=joinpath(expanduser(dir), "dftk_$(nm)_$stamp"),
                                   kspacing=parse(Float64, ks), ecutwfc=parse(Float64, ec), pseudo=ps)
            state[:results] = nothing
            v.status[] = "DFTK started in " * state[:job].dir
            follow_job!(v, controls)
        catch err
            err isa InterruptException && rethrow()
            v.status[] = "Could not start DFTK: " * sprint(showerror, err)
        end
    end
    Makie.on(stop.clicks) do _
        job = get(state, :job, nothing)
        job !== nothing && job_status(job) === :running || (v.status[] = "No DFTK job is running"; return)
        stop_job!(job)
        show_job!(v, controls)
        v.status[] = "DFTK job stopped"
    end
    Makie.on(open_result.clicks) do _
        r = get(state, :results, nothing)
        r === nothing && (v.status[] = "No results yet"; return)
        c = r.structure
        c = crystal_from_cartesian(c.lattice, c.species, cartesian_positions(c); pbc=c.pbc,
                                   title=c.title * (r.calculation == "relax" ? " (DFTK relaxed)" : " (DFTK)"))
        open_in_tab!(v, c)
        v.status[] = "Opened the $(r.calculation == "relax" ? "relaxed" : "calculated") structure (not saved yet)"
    end
    finish_dialog!(v, :dftk, "Run DFTK", f)
    return controls
end

busy(state) = (job = get(state, :job, nothing); job !== nothing && job_status(job) === :running)

"""Show the state of the tool's job: progress from the log, or the results when done."""
function show_job!(v::Viewer, t)
    job = get(t.state, :job, nothing)
    job === nothing && return :none
    st = job_status(job)
    secs = round(Int, time() - job.started)
    what = job.crystal === nothing ? "Installing DFTK" : "DFTK " * basename(job.script)
    tail = job_log(job; n=2)
    log = isempty(tail) ? "" : "\n" * join((first(strip(l), 60) for l in tail), "\n")
    if st === :running
        t.progress[] = "$what: running ($secs s)$log"
    elseif st === :done
        if job.crystal === nothing
            t.progress[] = "DFTK installed ($secs s)"
            t.envstate[] = "DFTK is installed there"
        else
            r = get(t.state, :results, nothing)
            r === nothing && (r = t.state[:results] = job_results(job))
            t.progress[] = "$what: done ($secs s), files in " * basename(job.dir) * "/"
            t.summary[] = summarize(r)
        end
    else
        t.progress[] = "$what: $(st === :stopped ? "stopped" : "failed, see " * basename(job.log))$log"
    end
    return st
end

"""Update the tool about once a second while the job runs (on drawn frames of the window)."""
function follow_job!(v::Viewer, t)
    stop = get(t.state, :follow, nothing)
    stop === nothing || stop()
    show_job!(v, t)
    last = Ref(0.0)
    t.state[:follow] = on_frames(v.figure.scene) do secs
        secs - last[] < 1 && return false
        last[] = secs
        st = show_job!(v, t)
        st === :running && return false
        job = t.state[:job]
        v.status[] = st === :done ? (job.crystal === nothing ? "DFTK is installed" :
                                     "DFTK finished: " * split(t.summary[], "\n")[2]) :
                     st === :failed ? "DFTK job failed; see " * job.log : v.status[]
        return true
    end
    return t
end
