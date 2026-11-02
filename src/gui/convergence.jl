# Plot of a DFTK convergence check (cutoff and k-grid series), its own window.

"""
    convergence_figure(cr::ConvergenceResults; tol=1.0, dir=nothing, name="convergence") -> NamedTuple

How the total energy converges: |E − E_best| in meV per atom (log scale)
against the cutoff (left) and against the k-spacing (right, coarse to
fine), with the tolerance `tol` dashed and the converged value
([`converged_settings`](@ref)) marked; the band gap of each step below when
there is one. With `dir`, a Save button writes `name.convergence.png` there.
Not displayed; see [`plot_convergence`](@ref).
"""
function convergence_figure(cr::ConvergenceResults; tol::Real=1.0, dir=nothing, name::AbstractString="convergence")
    fig = Makie.Figure(; size=(900, 560), backgroundcolor=:white)
    rec = converged_settings(cr; tol)
    heading = Makie.GridLayout(fig[0, 1:2])
    Makie.Label(heading[1, 1], "Convergence check"; font=:bold, fontsize=16)
    Makie.Label(heading[2, 1], "DFTK $(cr.version), $(cr.natoms) atoms, $(length(cr.steps)) of $(cr.planned) SCF runs; " *
                               "dashed: $(tol) meV/atom"; fontsize=13, color=:gray30)
    Makie.rowgap!(heading, 2)
    axes = Any[]
    hasgap = any(s -> s.gap !== nothing, cr.steps)
    for (col, series, xlabel, best) in ((1, "cutoff", "Cutoff (Ry)", rec.ecut), (2, "kpoints", "k-spacing (1/Å, incl. 2π)", rec.kspacing))
        s = convergence_series(cr, series)
        ax = Makie.Axis(fig[1, col]; xlabel=hasgap ? "" : xlabel, ylabel=col == 1 ? "|E − E_best| (meV/atom)" : "",
                        yscale=log10, xreversed=series == "kpoints", xgridvisible=false, yminorticksvisible=true,
                        title=series == "cutoff" ? "Cutoff" : "k-points")
        push!(axes, ax)
        if length(s.steps) < 2
            Makie.text!(ax, 0.5, 0.5; text=isempty(s.steps) ? "no runs" : "one k-grid only", space=:relative,
                        align=(:center, :center), color=:gray40)
            continue
        end
        err = [max(abs(e), 1e-4) for e in s.energy_error[1:end-1]]            # the last is the reference (0)
        x = s.parameter[1:end-1]
        Makie.scatterlines!(ax, x, err; color=DOS_LINE, markersize=9)
        Makie.hlines!(ax, [tol]; color=:gray30, linestyle=:dash)
        if best !== nothing && best != last(s.parameter)
            i = findfirst(==(best), x)
            Makie.scatter!(ax, [x[i]], [err[i]]; color=:firebrick, markersize=14, marker=:circle)
            series == "cutoff" && Makie.text!(ax, x[i], err[i]; text=@sprintf(" %g Ry", best), align=(:left, :bottom),
                                              color=:firebrick, fontsize=13)
        end
        if series == "kpoints"
            Makie.text!(ax, x, err; text=[join(st.kgrid, "×") for st in s.steps[1:end-1]], align=(:center, :top),
                        offset=(0, -8), fontsize=10, color=:gray40)
        end
        if hasgap
            axg = Makie.Axis(fig[2, col]; xlabel, ylabel=col == 1 ? "Gap (eV)" : "", xreversed=series == "kpoints",
                             xgridvisible=false)
            gx = [p for (p, g) in zip(s.parameter, s.gap) if g !== nothing]
            gy = [g for g in s.gap if g !== nothing]
            if !isempty(gx)
                Makie.scatterlines!(axg, gx, gy; color=DOS_DOWN_LINE, markersize=7)
                pad = max(0.05, 0.15 * (maximum(gy) - minimum(gy)))       # limits from the data (not displayed yet)
                Makie.ylims!(axg, max(0.0, minimum(gy) - pad), maximum(gy) + pad)
            end
            Makie.linkxaxes!(ax, axg)
            Makie.hidexdecorations!(ax; grid=false, ticks=false)
            push!(axes, axg)
        end
    end
    hasgap && Makie.rowsize!(fig.layout, 2, Makie.Relative(0.3))
    save = dir === nothing ? nothing : Makie.Button(fig[3, 1:2]; label="Save PNG to " * basename(dir) * "/")
    if save !== nothing
        Makie.on(save.clicks) do _
            Makie.save(joinpath(dir, name * ".convergence.png"), fig; update=false)
            save.label[] = "Saved $name.convergence.png"
        end
    end
    return (; figure=fig, axes, save, settings=rec)
end

"""
    plot_convergence(cr::ConvergenceResults; kw...) -> NamedTuple

Open the convergence plot ([`convergence_figure`](@ref)) in its own window.
"""
function plot_convergence(cr::ConvergenceResults; kw...)
    d = convergence_figure(cr; kw...)
    GLMakie.activate!()
    display(GLMakie.Screen(; title="DFTK convergence"), d.figure)
    return d
end
