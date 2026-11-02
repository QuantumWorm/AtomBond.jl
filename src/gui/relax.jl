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
    base = replace(c0.title, r" \(DFTK step .*\)$" => "")
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
    relax_figure(r::DFTKResults; on_step=nothing, dir=nothing, name="relax") -> NamedTuple

The steps of a relaxation ([`RelaxHistory`](@ref)): the energy relative to
the last step (meV) and the largest force on a free atom (eV/Å, log scale,
DFTK's 1 meV/Å target dashed) against the step. The *Step* slider marks a
step; `on_step(k, arrows)` is called when it or the *Force arrows* switch
changes (the Run DFTK tool shows that step in the main window). With
`on_bands` (called without arguments) a *Band structure / DOS* button is
added, for a run that computed them for the relaxed structure. With `dir`,
Save writes `name.relax.png` there.
"""
function relax_figure(r::DFTKResults; on_step=nothing, on_bands=nothing, dir=nothing, name::AbstractString="relax")
    h = r.history
    h === nothing && throw(ArgumentError("no relaxation steps in these results"))
    n = length(h)
    fig = Makie.Figure(; size=(760, 600), backgroundcolor=:white)
    heading = Makie.GridLayout(fig[0, 1])
    Makie.Label(heading[1, 1], "Relaxation of " * replace(h.structures[1].title, r" \(DFTK step .*\)$" => "");
                font=:bold, fontsize=16)
    Makie.Label(heading[2, 1], @sprintf("%d steps, energy change %.4f eV, %s", n - 1, last(h.energies) - first(h.energies),
                                        r.converged ? "converged" : "NOT converged"); fontsize=13, color=:gray30)
    Makie.rowgap!(heading, 2)
    steps = collect(0:n-1)
    axe = Makie.Axis(fig[1, 1]; ylabel="E − E_last (meV)", xgridvisible=false)
    Makie.scatterlines!(axe, steps, 1000 .* (h.energies .- last(h.energies)); color=DOS_LINE, markersize=7)
    axf = Makie.Axis(fig[2, 1]; xlabel="Step", ylabel="Largest force (eV/Å)", yscale=log10, xgridvisible=false,
                     xticks=collect(0:max(1, cld(n - 1, 10)):n-1))
    Makie.scatterlines!(axf, steps, max.(h.max_forces, 1e-5); color=DOS_DOWN_LINE, markersize=7)
    Makie.hlines!(axf, [1e-3]; color=:gray30, linestyle=:dash)
    Makie.linkxaxes!(axe, axf)
    Makie.hidexdecorations!(axe; grid=false, ticks=false)
    Makie.xlims!(axf, -0.5, n - 0.5)
    ev = 1000 .* (h.energies .- last(h.energies))
    Makie.ylims!(axe, min(minimum(ev), 0.0) - 0.05 * (maximum(ev) - minimum(ev) + 1), maximum(ev) * 1.05 + 0.5)
    Makie.ylims!(axf, max(minimum(h.max_forces), 1e-5) / 2, max(maximum(h.max_forces), 2e-3) * 2)
    controls = Makie.GridLayout(fig[3, 1])
    Makie.Label(controls[1, 1], "Step"; halign=:right)
    slider = Makie.Slider(controls[1, 2]; range=0:n-1, startvalue=n - 1)
    info = Makie.lift(k -> @sprintf("step %d: E − E_last = %.2f meV, largest force %.4f eV/Å", k,
                                    1000 * (h.energies[k+1] - last(h.energies)), h.max_forces[k+1]), slider.value)
    Makie.Label(controls[2, 1:4], info; halign=:left, fontsize=13)
    arrows = Makie.Toggle(controls[1, 3]; active=true)
    Makie.Label(controls[1, 4], "Force arrows"; halign=:left)
    Makie.colsize!(controls, 2, Makie.Relative(0.55))
    marker = Makie.lift(k -> [Float64(k)], slider.value)
    Makie.vlines!(axe, marker; color=:firebrick, linewidth=1)
    Makie.vlines!(axf, marker; color=:firebrick, linewidth=1)
    if on_step !== nothing
        Makie.onany(slider.value, arrows.active) do k, a
            on_step(k + 1, a)
        end
    end
    buttons = Makie.GridLayout(fig[4, 1])
    bands = on_bands === nothing ? nothing : Makie.Button(buttons[1, 1]; label="Band structure / DOS of the result")
    bands === nothing || Makie.on(_ -> on_bands(), bands.clicks)
    save = dir === nothing ? nothing : Makie.Button(buttons[1, 2]; label="Save PNG to " * basename(dir) * "/")
    if save !== nothing
        Makie.on(save.clicks) do _
            Makie.save(joinpath(dir, name * ".relax.png"), fig; update=false)
            save.label[] = "Saved $name.relax.png"
        end
    end
    return (; figure=fig, energy_axis=axe, force_axis=axf, step=slider, arrows, save, bands)
end

"""
    plot_relaxation(r::DFTKResults; kw...) -> NamedTuple

Open the relaxation steps ([`relax_figure`](@ref)) in their own window.
"""
function plot_relaxation(r::DFTKResults; kw...)
    d = relax_figure(r; kw...)
    GLMakie.activate!()
    display(GLMakie.Screen(; title="Relaxation — " * formula(r.structure)), d.figure)
    return d
end
