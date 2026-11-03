# Relaxation steps of a DFTK run (energy and forces per step, the structure of
# each step in a tab of the main window) and force arrows in the 3D view.

const FORCE_COLOR = RGBf(0.85, 0.12, 0.10)

"""
    show_forces!(viewer, forces; scale=nothing) -> scale

Draw the `forces` (eV/Å, one per atom of the current structure) as arrows
from the atoms in the 3D view, on every drawn copy of an atom. `scale` is
the arrow length in Å per eV/Å; by default the largest force gets a 2 Å
arrow. The arrows are part of the structure drawing, so they go away at the
next redraw (any change of the structure or of the view settings). Returns
the scale used.
"""
function show_forces!(v::Viewer, forces; scale=nothing)
    length(forces) == natoms(v.crystal) ||
        throw(ArgumentError("$(length(forces)) forces for $(natoms(v.crystal)) atoms"))
    fmax = maximum(norm, forces; init=0.0)
    s = scale === nothing ? (fmax > 0 ? 2.0 / fmax : 1.0) : Float64(scale)
    pts, dirs = Point3f[], Vec3f[]
    for a in v.display.atoms
        d = s .* forces[a.index]
        norm(d) < 0.02 && continue                       # too short to see
        push!(pts, Point3f(a.position...)); push!(dirs, Vec3f(d...))
    end
    if !isempty(pts)
        p = Makie.arrows3d!(v.lscene, pts, dirs; color=FORCE_COLOR, shaftradius=0.05, tipradius=0.12,
                            tiplength=0.3, markerscale=1)
        push!(v.plots, p)                                # removed with the structure plots at the next redraw
    end
    return s
end

"""Title marker of the tab that shows the relaxation steps."""
const RELAX_TAB_MARK = " (DFTK relaxation, step "

"""The number of the tab showing relaxation steps of `base`, or `nothing`."""
function relax_tab(v::Viewer, base::AbstractString)
    title(k) = k == v.active ? v.crystal.title : v.tabs[k].crystal.title
    return findfirst(k -> startswith(title(k), base * RELAX_TAB_MARK), 1:ntabs(v))
end

"""
Show step `k` of `h` in the relaxation tab of the main window (opened if
needed), with force arrows when `arrows`; the arrow scale is the same for
every step (the largest force of the whole relaxation gets 2 Å).
"""
function show_relax_step!(v::Viewer, h::RelaxHistory, k::Integer; arrows::Bool=true)
    c0 = h.structures[k]
    base = replace(c0.title, r" \(\w+ step \d+\)$" => "")
    c = Crystal(c0.lattice, c0.species, c0.frac; pbc=c0.pbc,
                fixed=[(x, x, x) for x in h.clamped], title=base * RELAX_TAB_MARK * "$(k - 1))")
    t = relax_tab(v, base)
    if t === nothing
        open_in_tab!(v, c)
    else
        switch_tab!(v, t)
        set_crystal!(v, c; reframe=false)
    end
    if arrows
        fmax = maximum(f -> maximum(norm, f; init=0.0), h.forces; init=0.0)
        show_forces!(v, h.forces[k]; scale=fmax > 0 ? 2.0 / fmax : 1.0)
    end
    return v
end

"""
    relax_figure(r::DFTResults; on_step=nothing, on_bands=nothing, dir=nothing, name="relax") -> NamedTuple
    relax_figure(history::Observable{RelaxHistory}; state=Observable("running"), kw...) -> NamedTuple

The steps of a relaxation ([`RelaxHistory`](@ref)): the energy relative to
the last step (meV) and the largest force on a free atom (eV/Å, log scale,
DFTK's 1 meV/Å target dashed) against the step; for a cell relaxation also
the volume change (%) and the pressure (GPa). The *Step* slider marks a
step; `on_step(k, arrows)` is called when it or the *Force arrows* switch
changes (the Run calculation tool shows that step in the main window). With
`on_bands` (called without arguments) a *Band structure / DOS* button is
added, for a run that computed them for the relaxed structure. With `dir`,
Save writes `name.relax.png` there.

Given an `Observable` history (a relaxation still running), the figure
follows it: new steps are added as the observable changes, and the slider
moves on to the newest step when it was on the last one; such automatic
moves call `on_new_step(k, arrows)` instead of `on_step`. `state` is shown in
the subtitle ("running", "converged", …).
"""
function relax_figure(r::DFTResults; kw...)
    r.history === nothing && throw(ArgumentError("no relaxation steps in these results"))
    return relax_figure(Makie.Observable(r.history); state=Makie.Observable(r.converged ? "converged" : "NOT converged"),
                        kw...)
end

function relax_figure(history::Makie.Observable{RelaxHistory}; state=Makie.Observable("running"), on_step=nothing,
                      on_new_step=nothing, on_bands=nothing, dir=nothing, name::AbstractString="relax")
    h0 = history[]
    fig = Makie.Figure(; size=(760, variable_cell(h0) ? 760 : 600), backgroundcolor=:white)
    heading = Makie.GridLayout(fig[0, 1]; tellwidth=false)
    start_title = replace(h0.structures[1].title, r" \(\w+ step \d+\)$" => "")
    Makie.Label(heading[1, 1], "Relaxation of " * (isempty(start_title) ? formula(h0.structures[1]) : start_title);
                font=:bold, fontsize=16)
    sub = Makie.lift((h, s) -> @sprintf("%s, energy change %.4f eV, %s", plural(length(h) - 1, "step"),
                                        last(h.energies) - first(h.energies), s), history, state)
    Makie.Label(heading[2, 1], sub; fontsize=13, color=:gray30)
    Makie.rowgap!(heading, 2)
    epts = Makie.lift(h -> [Makie.Point2f(k - 1, 1000 * (e - last(h.energies))) for (k, e) in enumerate(h.energies)], history)
    fpts = Makie.lift(h -> [Makie.Point2f(k - 1, max(f, 1e-5)) for (k, f) in enumerate(h.max_forces)], history)
    axe = Makie.Axis(fig[1, 1]; ylabel="E − E_last (meV)", xgridvisible=false)
    Makie.scatterlines!(axe, epts; color=DOS_LINE, markersize=7)
    axf = Makie.Axis(fig[2, 1]; xlabel="Step", ylabel="Largest force (eV/Å)", yscale=log10, xgridvisible=false)
    Makie.scatterlines!(axf, fpts; color=DOS_DOWN_LINE, markersize=7)
    Makie.hlines!(axf, [1e-3]; color=:gray30, linestyle=:dash)
    Makie.linkxaxes!(axe, axf)
    Makie.hidexdecorations!(axe; grid=false, ticks=false)
    cell = variable_cell(h0)
    axv = axp = nothing
    if cell                                   # volume change and pressure per step
        vpts = Makie.lift(h -> [Makie.Point2f(k - 1, cell_change(h, k).volume) for k in 1:length(h)], history)
        ppts = Makie.lift(h -> [Makie.Point2f(k - 1, p) for (k, p) in enumerate(pressure(h))], history)
        axv = Makie.Axis(fig[3, 1]; xlabel="Step", ylabel="Volume change (%)", xgridvisible=false)
        Makie.scatterlines!(axv, vpts; color=DOS_LINE, markersize=7)
        axp = Makie.Axis(fig[3, 1]; ylabel="Pressure (GPa)", yaxisposition=:right, ylabelcolor=DOS_DOWN_LINE,
                         yticklabelcolor=DOS_DOWN_LINE)
        Makie.hidexdecorations!(axp); Makie.hidespines!(axp); axp.ygridvisible = false
        Makie.scatterlines!(axp, ppts; color=DOS_DOWN_LINE, markersize=7, linestyle=:dash)
        Makie.hlines!(axp, [0.0]; color=:gray60, linestyle=:dot)
        Makie.linkxaxes!(axf, axv, axp)
        axf.xlabel = ""
        Makie.hidexdecorations!(axf; grid=false, ticks=false)
    end
    span(lo, hi; floor=0.0) = (d = max(hi - lo, floor); (lo - 0.08d, hi + 0.08d))
    function fit!(h)                          # limits from the data (autolimits needs a displayed figure)
        n = length(h)
        ev = 1000 .* (h.energies .- last(h.energies))
        ticks = collect(0:max(1, cld(n - 1, 10)):max(n - 1, 0))
        foreach(ax -> ax === nothing || (ax.xticks[] = ticks), (axf, axv, axp))
        Makie.xlims!(axf, -0.5, n - 0.5)
        if cell && variable_cell(h)
            vs = [cell_change(h, k).volume for k in 1:n]
            ps = pressure(h)
            Makie.xlims!(axv, -0.5, n - 0.5); Makie.xlims!(axp, -0.5, n - 0.5)
            Makie.ylims!(axv, span(min(minimum(vs), 0.0), max(maximum(vs), 0.0); floor=0.1)...)
            Makie.ylims!(axp, span(min(minimum(ps), 0.0), max(maximum(ps), 0.0); floor=0.1)...)
        end
        Makie.ylims!(axe, min(minimum(ev), 0.0) - 0.05 * (maximum(ev) - minimum(ev) + 1), maximum(ev) * 1.05 + 0.5)
        Makie.ylims!(axf, max(minimum(h.max_forces), 1e-5) / 2, max(maximum(h.max_forces), 2e-3) * 2)
    end
    fit!(h0)
    row = cell ? 4 : 3
    controls = Makie.GridLayout(fig[row, 1])
    Makie.Label(controls[1, 1], "Step"; halign=:right)
    n0 = length(h0)
    slider = Makie.Slider(controls[1, 2]; range=0:max(n0 - 1, 1), startvalue=n0 - 1)
    info = Makie.lift(slider.value, history) do k, h
        k = clamp(k, 0, length(h) - 1)
        t = @sprintf("step %d: E − E_last = %.2f meV, largest force %.4f eV/Å", k,
                     1000 * (h.energies[k+1] - last(h.energies)), h.max_forces[k+1])
        variable_cell(h) && (t *= @sprintf(", volume %+.2f %%, pressure %.3f GPa", cell_change(h, k + 1).volume,
                                           pressure(h.stresses[k+1])))
        t
    end
    Makie.Label(controls[2, 1:4], info; halign=:left, fontsize=13, tellwidth=false)
    arrows = Makie.Toggle(controls[1, 3]; active=true)
    Makie.Label(controls[1, 4], "Force arrows"; halign=:left)
    Makie.colsize!(controls, 2, Makie.Relative(0.6))
    controls.tellwidth = false
    marker = Makie.lift(k -> [Float64(k)], slider.value)
    Makie.vlines!(axe, marker; color=:firebrick, linewidth=1)
    Makie.vlines!(axf, marker; color=:firebrick, linewidth=1)
    cell && Makie.vlines!(axv, marker; color=:firebrick, linewidth=1)
    automatic = Ref(false)                    # the slider is being moved to a new step, not by the user
    last_n = Ref(n0)
    Makie.on(history) do h
        n = length(h)
        fit!(h)
        follow = slider.value[] >= last_n[] - 1
        slider.range[] = 0:max(n - 1, 1)
        if follow && n - 1 != slider.value[]
            automatic[] = true
            Makie.set_close_to!(slider, n - 1)
            automatic[] = false
            on_new_step === nothing || on_new_step(n, arrows.active[])
        end
        last_n[] = n
    end
    if on_step !== nothing
        Makie.onany(slider.value, arrows.active) do k, a
            automatic[] && return
            k + 1 <= length(history[]) && on_step(k + 1, a)
        end
    end
    # the buttons row only when there is a button (an empty row would leave a blank band)
    buttons = on_bands === nothing && dir === nothing ? nothing : Makie.GridLayout(fig[row+1, 1]; tellwidth=false)
    bands = on_bands === nothing ? nothing : Makie.Button(buttons[1, 1]; label="Band structure / DOS of the result")
    bands === nothing || Makie.on(_ -> on_bands(), bands.clicks)
    save = dir === nothing ? nothing : Makie.Button(buttons[1, 2]; label="Save PNG to " * basename(dir) * "/")
    if save !== nothing
        Makie.on(save.clicks) do _
            Makie.save(joinpath(dir, name * ".relax.png"), fig; update=false)
            save.label[] = "Saved $name.relax.png"
        end
    end
    return (; figure=fig, energy_axis=axe, force_axis=axf, volume_axis=axv, pressure_axis=axp, step=slider, arrows, save, bands, history, state)
end

"""
    plot_relaxation(r; kw...) -> NamedTuple

Open the relaxation steps ([`relax_figure`](@ref); `r` is a `DFTResults` or
an `Observable{RelaxHistory}`) in their own window.
"""
function plot_relaxation(r; kw...)
    d = relax_figure(r; kw...)
    GLMakie.activate!()
    display(GLMakie.Screen(; title="Relaxation — " * formula(d.history[].structures[1])), d.figure)
    return d
end
