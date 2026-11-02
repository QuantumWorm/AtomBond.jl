# Run calculation tool: start a calculation of the current structure with a DFT
# engine (DFTK by default; Quantum ESPRESSO, CASTEP, …), follow its progress and
# show the results when it is done (src/run/ does the work).

"""Tasks of the Run calculation tool (each engine offers some, see `engine_tasks`)."""
const DFTK_TASKS = [("SCF + band structure", :bands), ("SCF only", :scf), ("Relax atoms", :relax_atoms),
                    ("Relax atoms and cell", :relax_cell), ("Convergence check", :convergence)]

"""The `relax` setting of [`run_dftk`](@ref) for a task of the tool."""
dftk_relax_mode(task::Symbol) = task === :relax_atoms ? :atoms : task === :relax_cell ? :cell : :auto

"""The tasks of the tool `engine` can do, as menu options."""
engine_task_options(engine) = filter(t -> t[2] in engine_tasks(dft_engine(engine)), DFTK_TASKS)

"""The calculation the tool will set up for `c`, in words."""
function dftk_plan(c::Crystal, kspacing, task::Symbol=:scf; engine::Symbol=:dftk)
    if task === :convergence
        k = kspacing === nothing ? "" : " (k-grid " * join(kpoint_grid(c; kspacing), "×") * ")"
        nk = length(convergence_kspacings(c))
        return "Convergence check, $(plural(natoms(c), "atom")): 7 cutoffs$k,\n" *
               (nk > 1 ? "then $nk k-grids at the cutoff below" : "one k-point only: no k-grid series")
    end
    k = kspacing === nothing ? "" : ", k-grid " * join(kpoint_grid(c; kspacing), " × ")
    if task === :relax_atoms
        nfix = count(all, c.fixed)
        return "Relaxation of the atoms (cell fixed)\n" * plural(natoms(c), "atom") * k *
               (nfix > 0 ? "\n" * plural(nfix, "atom") * " kept fixed" : "")
    elseif task === :relax_cell
        warn = any(vacuum_axes(c)) || !all(c.pbc) ? "\n⚠ vacuum in the cell: it would shrink" :
               c.cell_constraints != FREE_CELL && c.cell_constraints != (0, 0, 0, 0, 0, 0) ?
                   (engine === :dftk ? "\nWhole cell relaxes (cell constraints not applied)" :
                    "\nCell: " * describe_cell_constraints(c.cell_constraints)) :
               engine === :dftk ? "\nStress needs a converged cutoff (use high)" : "\nStress needs a converged cutoff"
        return "Relaxation of the atoms and the cell\n" * plural(natoms(c), "atom") * k * warn
    end
    bands = task === :bands
    relax = has_constraints(c)
    what = relax ? "Relaxation (the structure has constraints)" : "SCF (energies, forces, levels)"
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

Run calculation tool: runs the current structure with a DFT engine in a
separate process. *Engine*: DFTK.jl (the default, in its own Julia
environment) or another registered engine ([`dft_engines`](@ref): Quantum
ESPRESSO, CASTEP, or one you add); for those the *Command* field shows
the command found ([`find_engine`](@ref)) and *Save settings* stores it,
with the pseudopotential folder, in the settings file
([`engine_settings_file`](@ref)). *Task*: an SCF with or without the band structure and dense DOS
([`run_dftk`](@ref); a relaxation when the structure has constraints), a
relaxation of the atoms or of the atoms and the cell, or a
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
When the environment has no DFTK yet, *Run* installs it first
([`install_dftk`](@ref), minutes, once) and then starts the calculation;
*Install DFTK* only installs it.
"""
function dftk_dialog(v::Viewer)
    f = new_dialog(v, :dftk, "Run calculation")
    f === nothing && return nothing
    heading!(f, "DFT calculation")
    plan = Makie.Observable("")
    Makie.Label(f.g[nextrow!(f), 1:2], plan; halign=:left, fontsize=13, color=:gray25, justification=:left)

    small_button(r, k, label) = with_click_feedback!(Makie.Button(f.g[r, k]; label, width=Makie.Relative(1),
                                                                  fontsize=PANEL_FONTSIZE, padding=(6, 6, 5, 5)))
    state = get!(v.dialogs, :dftk_state, Dict{Symbol,Any}())    # the job outlives the tool
    heading!(f, "Program")
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], "Engine"; halign=:left, fontsize=PANEL_FONTSIZE)
    engines = [(engine_name(e), engine_id(e)) for e in dft_engines()]
    chosen = get(state, :engine, :dftk)
    engine = Makie.Menu(f.g[r, 2]; options=engines, i_selected=something(findfirst(e -> e[2] === chosen, engines), 1),
                        fontsize=PANEL_FONTSIZE, width=Makie.Relative(1))
    eng() = engine.selection[]
    envlabel = Makie.Observable("Julia environment")
    env = field!(f, envlabel, default_dftk_env(), s -> !isempty(strip(s)))
    envstate = Makie.Observable("")
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], envstate; halign=:left, fontsize=12, color=:gray25, justification=:left)
    install = small_button(r, 2, "Install DFTK")

    heading!(f, "Calculation")
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], "Task"; halign=:left, fontsize=PANEL_FONTSIZE)
    task = Makie.Menu(f.g[r, 2]; options=engine_task_options(eng()), i_selected=1, fontsize=PANEL_FONTSIZE,
                      width=Makie.Relative(1))
    folder = field!(f, "Folder", pwd(), s -> !isempty(strip(s)))
    name = field!(f, "Name", formula(v.crystal), s -> occursin(r"^[\w.+-]+$", strip(s)))
    kspacing = field!(f, "k-spacing (1/Å, incl. 2π)", "0.25", s -> isnum(s) && parse(Float64, s) > 0)
    ecut = field!(f, "Cutoff (Ry)", "40", s -> isnum(s) && parse(Float64, s) > 0)
    hint = Makie.Observable("")
    r = nextrow!(f)
    Makie.Label(f.g[r, 1], hint; halign=:left, fontsize=12, color=:gray25, justification=:left)
    use_hint = small_button(r, 2, "Use recommended")
    pseudo = field!(f, "Pseudopotentials", engine_pseudo_default(dft_engine(eng())),
                    s -> eng() !== :dftk ? !isempty(strip(s)) : occursin("{El}", s) || occursin(r"^[\w.+-]+$", s))
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
    update_plan(_...) = (plan[] = (eng() === :dftk ? "" : engine_name(dft_engine(eng())) * ": ") *
                                  dftk_plan(v.crystal, num(kspacing), something(task.selection[], :scf); engine=eng()))
    Makie.on(update_plan, kspacing.displayed_string)
    Makie.on(update_plan, task.selection)
    Makie.on(update_plan, v.header)                  # the structure changed
    update_plan()
    function update_env(_...)
        if eng() === :dftk
            envstate[] = engine_status(DFTKEngine(); env=strip(shown_text(env)))
        else
            typed = strip(shown_text(env))
            f = find_engine(eng())
            envstate[] = isempty(typed) ? engine_status(eng()) :
                         f !== nothing && typed == join(f.command, " ") ? engine_status(eng()) :
                         (command_info(typed, "").found ? "Command typed in\n(Save settings keeps it)" :
                          "⚠ program not found:\ncheck the command")
        end
    end
    Makie.on(update_env, env.displayed_string)
    update_env()

    controls = (; plan, engine, env, envlabel, envstate, install, task, folder, name, kspacing, ecut, hint, use_hint, pseudo, run, stop,
                progress, summary, open_result, show_dos, state)

    # Cutoff hints: from the DFTK environment, fetched in the background (a helper process the first time).
    hints = Makie.Observable(Dict{Symbol,Any}())
    update_hint(_...) = (hint[] = eng() === :dftk ? cutoff_hint_text(v.crystal, hints[], num(ecut)) :
                                  "No cutoff hints for\n" * engine_name(dft_engine(eng())) * ": check the pseudos")
    function fetch_hints(_...)
        eng() === :dftk || return
        ccall(:jl_generating_output, Cint, ()) == 1 && return     # no helper process while precompiling
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
    # Another engine: its command (or DFTK's environment), tasks and pseudopotentials.
    function engine_changed(id)
        state[:engine] = id
        if id === :dftk
            envlabel[] = "Julia environment"
            settext!(env, default_dftk_env())
            install.label[] = "Install DFTK"
        else
            envlabel[] = "Command"
            fnd = find_engine(id)
            settext!(env, fnd === nothing ? "" : join(fnd.command, " "))
            install.label[] = "Save settings"
        end
        settext!(pseudo, engine_pseudo_default(dft_engine(id)))
        current = task.selection[]
        task.options[] = engine_task_options(id)
        k = findfirst(o -> o[2] === current, task.options[])
        task.i_selected[] = something(k, 1)
        update_env(); update_hint(); update_plan(); fetch_hints()
    end
    Makie.on(engine_changed, engine.selection)
    eng() === :dftk || engine_changed(eng())
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
        if eng() !== :dftk                         # Save settings: the command and the pseudopotentials
            cmd = form_value(v, env, "the command that runs $(engine_name(dft_engine(eng())))")
            cmd === nothing && return
            set_engine_setting!(eng(), "command", cmd)
            ps = strip(shown_text(pseudo))
            isempty(ps) || set_engine_setting!(eng(), "pseudo_dir", ps)
            update_env()
            v.status[] = "Saved the $(engine_name(dft_engine(eng()))) settings in " * engine_settings_file()
            return
        end
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
        id = eng()
        vals = (form_value(v, env, id === :dftk ? "an environment folder" : "the command that runs the program"),
                form_value(v, folder, "a folder"),
                form_value(v, name, "a name"), form_value(v, kspacing, "a k-point spacing"),
                form_value(v, ecut, "a cutoff"), form_value(v, pseudo, "pseudopotentials"))
        any(isnothing, vals) && return
        busy(state) && (v.status[] = "A calculation is still running: stop it first"; return)
        e, dir, nm, ks, ec, ps = vals
        c, t = v.crystal, task.selection[]
        if id !== :dftk
            ename = engine_name(dft_engine(id))
            try
                stamp = Dates.format(Dates.now(), "yyyymmdd-HHMMSS")
                state[:job] = run_dft(c; engine=id, command=e, name=nm, dir=joinpath(expanduser(dir), "$(id)_$(nm)_$stamp"),
                                      kspacing=parse(Float64, ks), ecutwfc=parse(Float64, ec), pseudo=ps,
                                      relax=dftk_relax_mode(t))
                state[:results] = nothing
                delete!(state, :live); delete!(state, :live_state)
                summary[] = ""
                v.status[] = "$ename started in " * state[:job].dir
                follow_job!(v, controls)
            catch err
                err isa InterruptException && rethrow()
                v.status[] = "Could not start $ename: " * sprint(showerror, err)
            end
            return
        end
        typed_ecut = strip(shown_text(ecut)) != auto_ecut[]
        function start()
            try
                stamp = Dates.format(Dates.now(), "yyyymmdd-HHMMSS")
                ec = typed_ecut ? ec : something(num(ecut), parse(Float64, ec))     # the recommended one, once known
                state[:job] = t === :convergence ?
                    run_dftk_convergence(c; env=e, name=nm, dir=joinpath(expanduser(dir), "dftk_$(nm)_conv_$stamp"),
                                         kspacing=parse(Float64, ks), ecutwfc=Float64(ec), pseudo=ps) :
                    run_dftk(c; env=e, name=nm, dir=joinpath(expanduser(dir), "dftk_$(nm)_$stamp"),
                             kspacing=parse(Float64, ks), ecutwfc=Float64(ec), pseudo=ps, bands=t === :bands,
                             relax=dftk_relax_mode(t))
                state[:results] = nothing
                delete!(state, :live); delete!(state, :live_state)
                summary[] = ""
                v.status[] = "DFTK started in " * state[:job].dir
                follow_job!(v, controls)
            catch err
                err isa InterruptException && rethrow()
                v.status[] = "Could not start DFTK: " * sprint(showerror, err)
            end
        end
        if dftk_env_ready(e)
            start()
        else                       # DFTK first (once), then the calculation; the cutoff hints come with it
            state[:job] = install_dftk(e)
            state[:results] = nothing
            state[:after_install] = function ()
                @async begin
                    try
                        hints[] = Dict{Symbol,Any}(dftk_cutoff_hints(species_list(c); family=ps, env=e))
                    catch err
                        err isa InterruptException && rethrow()
                    end
                    start()
                end
            end
            summary[] = "DFTK is not installed in $(basename(e)) yet:\ninstalling it first (a few minutes, once),\nthen the calculation starts by itself"
            v.status[] = "Installing DFTK into $e, then running $nm"
            follow_job!(v, controls)
        end
    end
    Makie.on(stop.clicks) do _
        job = get(state, :job, nothing)
        job !== nothing && job_status(job) === :running || (v.status[] = "No calculation is running"; return)
        stop_job!(job)
        show_job!(v, controls)
        v.status[] = "Calculation stopped"
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
                                   title=c.title * " (" * (isempty(r.program) ? "DFTK" : r.program) *
                                         (r.calculation == "relax" ? " relaxed)" : ")"))
        open_in_tab!(v, c)
        v.status[] = "Opened the $(r.calculation == "relax" ? "relaxed" : "calculated") structure (not saved yet)"
    end
    Makie.on(show_dos.clicks) do _
        r = get(state, :results, nothing)
        job = get(state, :job, nothing)
        if r === nothing && job !== nothing && is_relax_job(job) && job_status(job) === :running
            live = get(state, :live, nothing)
            live === nothing && (v.status[] = "No relaxation step yet: wait for the first one"; return)
            base = replace(splitext(basename(job.script))[1], r"\.conv$" => "")
            state[:dos] = (v.screen === nothing ? relax_figure : plot_relaxation)(live; state=state[:live_state],
                dir=job.dir, name=base,
                on_step=(k, arrows) -> show_relax_step!(v, live[], k; arrows),
                on_new_step=(k, arrows) -> relax_tab_active(v, live[]) && show_relax_step!(v, live[], k; arrows))
            v.status[] = "Relaxation in progress" * (v.screen === nothing ? "" : ", shown in a new window") *
                         ": it follows the new steps"
            return
        end
        r === nothing && (v.status[] = "No results yet"; return)
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
    finish_dialog!(v, :dftk, "Run calculation", f)
    return controls
end

"""Open the band structure and DOS (or the DOS alone) of `r` in its own window (a figure only without a window)."""
function show_levels!(v::Viewer, r::DFTResults, dir, base)
    state = v.dialogs[:dftk_state]
    withbands = r.path !== nothing
    show = withbands ? (v.screen === nothing ? bands_figure : plot_bands) : (v.screen === nothing ? dos_figure : plot_dos)
    state[:levels] = show(r; dir, name=base)
    v.status[] = (withbands ? "Band structure and DOS of " : "DOS of ") * formula(r.structure) *
                 (v.screen === nothing ? "" : " opened in a new window")
end

busy(state) = (job = get(state, :job, nothing); job !== nothing && job_status(job) === :running)

"""Whether `job` is a relaxation (not a convergence check)."""
is_relax_job(job::DFTJob) = job.crystal !== nothing && job.relax && !is_convergence_job(job)

"""Whether the main window shows the relaxation tab of `h`."""
function relax_tab_active(v::Viewer, h::RelaxHistory)
    base = replace(h.structures[1].title, r" \(\w+ step \d+\)$" => "")
    return relax_tab(v, base) == v.active
end

"""Whether `job` is a convergence check."""
is_convergence_job(job::DFTJob) = endswith(job.results, ".convergence.toml")

"""Show the state of the tool's job: progress from the log, or the results when done."""
function show_job!(v::Viewer, t)
    job = get(t.state, :job, nothing)
    job === nothing && return :none
    st = job_status(job)
    secs = round(Int, time() - job.started)
    what = job.crystal === nothing ? "Installing DFTK" : engine_name(dft_engine(job.engine)) * " " * basename(job.script)
    tail = job_log(job; n=2)
    log = isempty(tail) ? "" : "\n" * join((first(strip(l), 60) for l in tail), "\n")
    conv = job.crystal !== nothing && is_convergence_job(job)
    relaxed = !conv && get(t.state, :results, nothing) isa DFTResults && t.state[:results].history !== nothing
    t.open_result.label[] = conv ? "Use converged values" : "Open structure"
    t.show_dos.label[] = conv ? "Show plot" : relaxed ? "Show relaxation" : "Show DOS / bands"
    if st === :running
        t.progress[] = "$what: running ($secs s)$log"
        if conv                                                     # partial results as they come
            cr = try job_convergence(job) catch; nothing end
            cr === nothing || (t.summary[] = summarize(cr))
        elseif is_relax_job(job)                                    # relaxation steps as they come
            h = job_history(job)
            if h !== nothing
                live = get!(() -> Makie.Observable(h), t.state, :live)
                get!(() -> Makie.Observable("running"), t.state, :live_state)
                length(h) != length(live[]) && (live[] = h)
                t.summary[] = "Relaxation running\n" * relax_progress(h)
                t.show_dos.label[] = "Show relaxation"
            end
        end
    elseif st === :done
        if job.crystal === nothing
            t.progress[] = "DFTK installed ($secs s)"
            t.envstate[] = "DFTK is installed there"
        else
            r = get(t.state, :results, nothing)
            r === nothing && (r = t.state[:results] = conv ? job_convergence(job) : job_results(job))
            if r isa DFTResults && r.history !== nothing && haskey(t.state, :live)   # a live window: the final steps
                length(r.history) != length(t.state[:live][]) && (t.state[:live][] = r.history)
                t.state[:live_state][] = r.converged ? "converged" : "NOT converged"
            end
            t.progress[] = "$what: done ($secs s), files in " * basename(job.dir) * "/"
            t.summary[] = summarize(r)
            r isa DFTResults && r.history !== nothing && (t.show_dos.label[] = "Show relaxation")
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
        after = pop!(t.state, :after_install, nothing)              # a calculation waiting for the installation
        if job.crystal === nothing && after !== nothing
            st === :done ? after() : (t.summary[] = "")
        end
        lines = split(t.summary[], "\n")
        ename = job.crystal === nothing ? "DFTK" : engine_name(dft_engine(job.engine))
        v.status[] = st === :done ? (job.crystal === nothing ? "DFTK is installed" :
                                     "$ename finished: " * lines[min(2, end)]) :
                     st === :failed ? "$ename failed; see " * job.log : v.status[]
        return true
    end
    return t
end
