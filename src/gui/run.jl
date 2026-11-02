# Run DFTK tool: start a DFTK calculation of the current structure, follow its
# progress and show the results when it is done (src/run/dftk.jl does the work).

"""Tasks of the Run DFTK tool."""
const DFTK_TASKS = [("SCF + band structure", :bands), ("SCF only", :scf), ("Convergence check", :convergence)]

"""The calculation the tool will set up for `c`, in words."""
function dftk_plan(c::Crystal, kspacing, task::Symbol=:scf)
    if task === :convergence
        k = kspacing === nothing ? "" : " (k-grid " * join(kpoint_grid(c; kspacing), "×") * ")"
        nk = length(convergence_kspacings(c))
        return "Convergence check, $(plural(natoms(c), "atom")): 7 cutoffs$k,\n" *
               (nk > 1 ? "then $nk k-grids at the cutoff below" : "one k-point only: no k-grid series")
    end
    bands = task === :bands
    relax = has_constraints(c)
    what = relax ? "Relaxation (the structure has constraints)" : "SCF (energies, forces, levels)"
    k = kspacing === nothing ? "" : ", k-grid " * join(kpoint_grid(c; kspacing), " × ")
    extra = !bands ? "" :
            (relax ? "\nThen an SCF of the result with " : "\nThen ") * "DOS on a " *
            (kspacing === nothing ? "denser grid" : join(kpoint_grid(c; kspacing=kspacing / 2), " × ")) *
            (any(vacuum_axes(c)) || !all(c.pbc) ? " (no band path: vacuum)" : " and bands")
    return what * "\n" * plural(natoms(c), "atom") * k * extra
end

"""The cutoff hint line for `c`: the recommended cutoffs, and a warning below the low one."""
function cutoff_hint_text(c::Crystal, hints, cutoff)
    isempty(hints) && return "No recommendation\n(needs DFTK, a family)"
    els = species_list(c)
    all(e -> haskey(hints, e), els) || return "No recommendation for\n" * join(filter(e -> !haskey(hints, e), els), ", ")
    lvl(k) = 2 * maximum(getproperty(hints[e], k) for e in els)
    t = @sprintf("Recommended %g Ry\n(low %g, high %g)", lvl(:normal), lvl(:low), lvl(:high))
    cutoff !== nothing && cutoff < lvl(:low) && (t = @sprintf("⚠ below the low %g Ry\nRecommended %g Ry", lvl(:low), lvl(:normal)))
    return t
end

"""
    dftk_dialog(viewer) -> NamedTuple of the tool's controls

Run DFTK tool: runs the current structure in DFTK.jl in a separate Julia
process. *Task*: an SCF with or without the band structure and dense DOS
([`run_dftk`](@ref); a relaxation when the structure has constraints), or a
convergence check of the cutoff and k-grid ([`run_dftk_convergence`](@ref)).
The cutoff line shows what the pseudopotentials recommend for the elements
of the structure ([`recommended_ecut`](@ref)); the cutoff starts at that
value once it is known (unless one was typed in), and *Use recommended*
puts it back. While a job runs the tool shows the end of its log; when it is done, a
summary. For an SCF, *Open structure* opens the final (relaxed) structure in
a new tab and *Show DOS / bands* the DOS or band structure in its own
window; for a relaxation, *Show relaxation* shows energy and forces per step
([`relax_figure`](@ref)) and each step, with force arrows, in a tab; for a convergence check, *Use converged values* puts the converged
cutoff and k-spacing in the form and *Show plot* plots the check.
*Install DFTK* creates the DFTK environment the first time
([`install_dftk`](@ref)).
"""
function dftk_dialog(v::Viewer)
    f = new_dialog(v, :dftk, "Run DFTK")
    f === nothing && return nothing
    heading!(f, "DFTK.jl calculation")
    plan = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], plan; halign=:left, fontsize=13, color=:gray25, justification=:left)

    small_button(r, k, label) = with_click_feedback!(Makie.Button(f.g[r, k]; label, width=Makie.Relative(1),
                                                                  fontsize=PANEL_FONTSIZE, padding=(6, 6, 5, 5)))
    heading!(f, "Environment")
    env = field!(f, "Julia environment", default_dftk_env(), s -> !isempty(strip(s)))
    envstate = Makie.Observable("")
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], envstate; halign=:left, fontsize=12, color=:gray25, justification=:left)
    install = small_button(r, 2, "Install DFTK")

    heading!(f, "Calculation")
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], "Task"; halign=:left, fontsize=PANEL_FONTSIZE)
    task = Makie.Menu(f.g[r, 2]; options=DFTK_TASKS, i_selected=1, fontsize=PANEL_FONTSIZE, width=Makie.Relative(1))
    folder = field!(f, "Folder", pwd(), s -> !isempty(strip(s)))
    name = field!(f, "Name", formula(v.crystal), s -> occursin(r"^[\w.+-]+$", strip(s)))
    kspacing = field!(f, "k-spacing (1/Å, incl. 2π)", "0.25", s -> isnum(s) && parse(Float64, s) > 0)
    ecut = field!(f, "Cutoff (Ry)", "40", s -> isnum(s) && parse(Float64, s) > 0)
    hint = Makie.Observable("")
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], hint; halign=:left, fontsize=12, color=:gray25, justification=:left)
    use_hint = small_button(r, 2, "Use recommended")
    pseudo = field!(f, "Pseudopotentials", DFTK_PSEUDO_FAMILY, s -> occursin("{El}", s) || occursin(r"^[\w.+-]+$", s))
    r = nextrow!(f)
    run, stop = small_button(r, 1, "Run"), small_button(r, 2, "Stop")

    heading!(f, "Progress and results")
    progress = Makie.Observable("Not started")
    Makie.Label(f.g[nextrow!(f), 1:2], progress; halign=:left, fontsize=12, color=:gray25, justification=:left)
    summary = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], summary; halign=:left, fontsize=13, justification=:left)
    r = nextrow!(f)
    open_result, show_dos = small_button(r, 1, "Open structure"), small_button(r, 2, "Show DOS / bands")

    num(t) = (x = tryparse(Float64, strip(shown_text(t))); x === nothing || x <= 0 ? nothing : x)
    update_plan(_...) = (plan[] = dftk_plan(v.crystal, num(kspacing), task.selection[]))
    Makie.on(update_plan, kspacing.displayed_string)
    Makie.on(update_plan, task.selection)
    Makie.on(update_plan, v.header)                  # the structure changed
    update_plan()
    update_env(_...) = (envstate[] = dftk_env_ready(strip(shown_text(env))) ? "DFTK is installed there" :
                                     "No DFTK there yet:\ninstall it (minutes, once)")
    Makie.on(update_env, env.displayed_string)
    update_env()

    state = get!(v.dialogs, :dftk_state, Dict{Symbol,Any}())    # the job outlives the tool
    controls = (; plan, env, envstate, install, task, folder, name, kspacing, ecut, hint, use_hint, pseudo, run, stop,
                progress, summary, open_result, show_dos, state)

    # Cutoff hints: from the DFTK environment, fetched in the background (a helper process the first time).
    hints = Makie.Observable(Dict{Symbol,Any}())
    update_hint(_...) = (hint[] = cutoff_hint_text(v.crystal, hints[], num(ecut)))
    function fetch_hints(_...)
        e, ps = strip(shown_text(env)), strip(shown_text(pseudo))
        dftk_env_ready(e) || (hints[] = Dict{Symbol,Any}(); return)
        @async try
            hints[] = Dict{Symbol,Any}(dftk_cutoff_hints(species_list(v.crystal); family=ps, env=e))
        catch err
            err isa InterruptException && rethrow()
        end
    end
    Makie.on(update_hint, hints)
    Makie.on(update_hint, ecut.displayed_string)
    # The cutoff starts at the recommended value as soon as it is known (unless one was typed in).
    auto_ecut = Ref(strip(shown_text(ecut)))
    Makie.on(hints) do h
        els = species_list(v.crystal)
        (isempty(els) || !all(e -> haskey(h, e), els)) && return
        rec = @sprintf("%g", 2 * maximum(h[e].normal for e in els))
        strip(shown_text(ecut)) == auto_ecut[] || return
        settext!(ecut, rec)
        auto_ecut[] = rec
    end
    foreach(o -> Makie.on(fetch_hints, o), (env.stored_string, pseudo.stored_string, v.header))
    update_hint()
    fetch_hints()
    state[:hints] = hints
    show_job!(v, controls)

    Makie.on(use_hint.clicks) do _
        rec = let h = hints[], els = species_list(v.crystal)
            all(e -> haskey(h, e), els) ? 2 * maximum(h[e].normal for e in els) : nothing
        end
        rec === nothing && (v.status[] = "No recommendation for these elements and pseudopotentials"; return)
        settext!(ecut, @sprintf("%g", rec))
        v.status[] = @sprintf("Cutoff set to the recommended %g Ry", rec)
    end
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
            t = task.selection[]
            state[:job] = t === :convergence ?
                run_dftk_convergence(v.crystal; env=e, name=nm, dir=joinpath(expanduser(dir), "dftk_$(nm)_conv_$stamp"),
                                     kspacing=parse(Float64, ks), ecutwfc=parse(Float64, ec), pseudo=ps) :
                run_dftk(v.crystal; env=e, name=nm, dir=joinpath(expanduser(dir), "dftk_$(nm)_$stamp"),
                         kspacing=parse(Float64, ks), ecutwfc=parse(Float64, ec), pseudo=ps, bands=t === :bands)
            state[:results] = nothing
            summary[] = ""
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
        if r isa ConvergenceResults
            rec = converged_settings(r)
            rec.ecut === nothing && rec.kspacing === nothing &&
                (v.status[] = "Nothing converged within the ranges of the check"; return)
            rec.ecut === nothing || settext!(ecut, @sprintf("%g", rec.ecut))
            rec.kspacing === nothing || settext!(kspacing, @sprintf("%g", rec.kspacing))
            v.status[] = "Converged values in the form: " *
                         join(filter(!isempty, [rec.ecut === nothing ? "" : @sprintf("cutoff %g Ry", rec.ecut),
                                                rec.kspacing === nothing ? "" : @sprintf("k-spacing %g 1/Å", rec.kspacing)]), ", ")
            return
        end
        c = r.structure
        c = crystal_from_cartesian(c.lattice, c.species, cartesian_positions(c); pbc=c.pbc,
                                   title=c.title * (r.calculation == "relax" ? " (DFTK relaxed)" : " (DFTK)"))
        open_in_tab!(v, c)
        v.status[] = "Opened the $(r.calculation == "relax" ? "relaxed" : "calculated") structure (not saved yet)"
    end
    Makie.on(show_dos.clicks) do _
        r = get(state, :results, nothing)
        r === nothing && (v.status[] = "No results yet"; return)
        job = state[:job]
        base = replace(splitext(basename(job.script))[1], r"\.conv$" => "")
        try
            if r isa ConvergenceResults
                state[:dos] = (v.screen === nothing ? convergence_figure : plot_convergence)(r; dir=job.dir, name=base)
                v.status[] = "Convergence plot" * (v.screen === nothing ? "" : " opened in a new window")
                return
            end
            if r.history !== nothing                                   # a relaxation: its steps, and bands if computed
                on_bands = isempty(r.eigenvalues) ? nothing : () -> show_levels!(v, r, job.dir, base)
                state[:dos] = (v.screen === nothing ? relax_figure : plot_relaxation)(r; dir=job.dir, name=base, on_bands,
                    on_step=(k, arrows) -> show_relax_step!(v, r.history, k; arrows))
                v.status[] = "Relaxation steps" * (v.screen === nothing ? "" : " opened in a new window") *
                             ": the Step slider shows each step in a tab"
                return
            end
            isempty(r.eigenvalues) && (v.status[] = "A relaxation keeps no levels: run an SCF of the result for the DOS"; return)
            show_levels!(v, r, job.dir, base)
        catch err
            err isa InterruptException && rethrow()
            v.status[] = "Could not show the plot: " * sprint(showerror, err)
        end
    end
    finish_dialog!(v, :dftk, "Run DFTK", f)
    return controls
end

"""Open the band structure and DOS (or the DOS alone) of `r` in its own window (a figure only without a window)."""
function show_levels!(v::Viewer, r::DFTKResults, dir, base)
    state = v.dialogs[:dftk_state]
    withbands = r.path !== nothing
    show = withbands ? (v.screen === nothing ? bands_figure : plot_bands) : (v.screen === nothing ? dos_figure : plot_dos)
    state[:levels] = show(r; dir, name=base)
    v.status[] = (withbands ? "Band structure and DOS of " : "DOS of ") * formula(r.structure) *
                 (v.screen === nothing ? "" : " opened in a new window")
end

busy(state) = (job = get(state, :job, nothing); job !== nothing && job_status(job) === :running)

"""Whether `job` is a convergence check."""
is_convergence_job(job::DFTKJob) = endswith(job.results, ".convergence.toml")

"""Show the state of the tool's job: progress from the log, or the results when done."""
function show_job!(v::Viewer, t)
    job = get(t.state, :job, nothing)
    job === nothing && return :none
    st = job_status(job)
    secs = round(Int, time() - job.started)
    what = job.crystal === nothing ? "Installing DFTK" : "DFTK " * basename(job.script)
    tail = job_log(job; n=2)
    log = isempty(tail) ? "" : "\n" * join((first(strip(l), 60) for l in tail), "\n")
    conv = job.crystal !== nothing && is_convergence_job(job)
    relaxed = !conv && get(t.state, :results, nothing) isa DFTKResults && t.state[:results].history !== nothing
    t.open_result.label[] = conv ? "Use converged values" : "Open structure"
    t.show_dos.label[] = conv ? "Show plot" : relaxed ? "Show relaxation" : "Show DOS / bands"
    if st === :running
        t.progress[] = "$what: running ($secs s)$log"
        if conv                                                     # partial results as they come
            cr = try job_convergence(job) catch; nothing end
            cr === nothing || (t.summary[] = summarize(cr))
        end
    elseif st === :done
        if job.crystal === nothing
            t.progress[] = "DFTK installed ($secs s)"
            t.envstate[] = "DFTK is installed there"
        else
            r = get(t.state, :results, nothing)
            r === nothing && (r = t.state[:results] = conv ? job_convergence(job) : job_results(job))
            t.progress[] = "$what: done ($secs s), files in " * basename(job.dir) * "/"
            t.summary[] = summarize(r)
            r isa DFTKResults && r.history !== nothing && (t.show_dos.label[] = "Show relaxation")
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
        lines = split(t.summary[], "\n")
        v.status[] = st === :done ? (job.crystal === nothing ? "DFTK is installed" :
                                     "DFTK finished: " * lines[min(2, end)]) :
                     st === :failed ? "DFTK job failed; see " * job.log : v.status[]
        return true
    end
    return t
end
