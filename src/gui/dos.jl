# Density-of-states plot of a finished DFTK run (its own window).

const DOS_FILL = RGBAf(0.20, 0.45, 0.75, 0.30)    # filled states
const DOS_LINE = RGBf(0.12, 0.30, 0.55)
const DOS_DOWN_FILL = RGBAf(0.80, 0.35, 0.20, 0.30)
const DOS_DOWN_LINE = RGBf(0.60, 0.22, 0.10)

"""
    dos_figure(r::DFTKResults; width=0.1, title="", dir=nothing, name="dos") -> (; figure, axis, width, save, ...)

A figure with the density of states of `r` ([`dos`](@ref)) against the
energy relative to the Fermi level: the filled states (below the Fermi
level) shaded, the Fermi level dashed, spin down mirrored below the axis
for a spin-polarized run. A slider sets the Gaussian broadening (eV) and the
plot follows it. Deep levels (semicore states, separated from the valence
levels by an empty range of more than `deep_gap` eV, 5 by default) are left
out of the view unless *Show deep levels* is on; the CSV always has them. With `dir`, a Save button writes `name.dos.png`
and `name.dos.csv` there. The figure is not displayed; see [`plot_dos`](@ref).
"""
function dos_figure(r::DFTKResults; width::Real=0.1, title::AbstractString="", dir=nothing, name::AbstractString="dos",
                    deep_gap::Real=5)
    isempty(r.eigenvalues) && throw(ArgumentError("no levels in these results (a relaxation keeps none)"))
    fig = Makie.Figure(; size=(760, 520), backgroundcolor=:white)
    ef = something(r.fermi, 0.0)
    head = !isempty(title) ? title : !isempty(r.structure.title) ? r.structure.title : formula(r.structure)
    nk = length(r.eigenvalues) ÷ r.n_spin
    sub = "DFTK $(r.version), $(nk) k-point$(nk == 1 ? "" : "s")" * (r.n_spin == 2 ? ", spin-polarized" : "") *
          (r.fermi === nothing ? "" : ", energies relative to the Fermi level")
    ax = Makie.Axis(fig[1, 1:3]; title=head, subtitle=sub, xlabel=r.fermi === nothing ? "E (eV)" : "E − E_F (eV)",
                    ylabel="DOS (states/eV per cell)", xgridvisible=false, ygridvisible=false)
    slider = Makie.Slider(fig[2, 2]; range=0.02:0.01:0.5, startvalue=clamp(width, 0.02, 0.5))
    Makie.Label(fig[2, 1], "Broadening"; halign=:right)
    wlabel = Makie.lift(w -> @sprintf("%.2f eV", w), slider.value)
    Makie.Label(fig[2, 3], wlabel; halign=:left, width=70)
    valence_bottom = deep_levels_edge(r; gap=deep_gap)        # nothing when there are no deep levels
    showdeep = Makie.Observable(false)
    if valence_bottom !== nothing
        gl = Makie.GridLayout(fig[3, 1:3]; tellwidth=false, halign=:left)
        tg = Makie.Toggle(gl[1, 1]; active=false)
        Makie.Label(gl[1, 2], @sprintf("Show deep levels (below %.0f eV)", valence_bottom - ef); halign=:left)
        Makie.connect!(showdeep, tg.active)
    end
    row = valence_bottom === nothing ? 3 : 4
    save = dir === nothing ? nothing : Makie.Button(fig[row, 1:3]; label="Save PNG and CSV to " * basename(dir) * "/")
    data = Makie.lift(w -> dos(r; width=w), slider.value)
    x = Makie.lift(d -> d.energies .- ef, data)
    for σ in 1:r.n_spin
        sgn = σ == 1 ? 1 : -1
        y = Makie.lift(d -> sgn .* d.dos[:, σ], data)
        filled = Makie.lift((xs, ys) -> [xi <= 0 ? yi : 0.0 for (xi, yi) in zip(xs, ys)], x, y)
        Makie.band!(ax, x, Makie.lift(v -> zero(v), filled), filled; color=σ == 1 ? DOS_FILL : DOS_DOWN_FILL)
        Makie.lines!(ax, x, y; color=σ == 1 ? DOS_LINE : DOS_DOWN_LINE, linewidth=1.5,
                     label=r.n_spin == 1 ? "DOS" : (σ == 1 ? "spin up" : "spin down"))
    end
    r.n_spin == 2 && Makie.hlines!(ax, [0.0]; color=:gray50, linewidth=0.8)
    r.fermi === nothing || Makie.vlines!(ax, [0.0]; color=:gray30, linestyle=:dash, linewidth=1)
    r.n_spin == 2 && Makie.axislegend(ax; position=:rt, framevisible=false)
    function fit!(_...)                   # limits from the data in view (autolimits needs a displayed figure)
        d = data[]
        xs = d.energies .- ef
        lo = showdeep[] || valence_bottom === nothing ? first(xs) : valence_bottom - ef - 4 * slider.value[]
        top = maximum(d.dos[xs .>= lo, :]) * 1.08
        Makie.limits!(ax, lo, last(xs), r.n_spin == 2 ? -top : 0.0, top)
    end
    Makie.on(fit!, data)
    Makie.on(fit!, showdeep)
    fit!()
    if save !== nothing
        Makie.on(save.clicks) do _
            png, csv = joinpath(dir, name * ".dos.png"), joinpath(dir, name * ".dos.csv")
            write_dos(csv, r; width=slider.value[])
            Makie.save(png, fig; update=false)
            save.label[] = "Saved " * basename(png) * " and " * basename(csv)
        end
    end
    return (; figure=fig, axis=ax, width=slider, save, data, showdeep)
end

"""
Lowest valence level (eV) when deep (semicore) levels lie below it, separated
by an empty range wider than `gap` eV below the Fermi level; otherwise `nothing`.
"""
function deep_levels_edge(r::DFTKResults; gap::Real=5)
    ef = something(r.fermi, Inf)
    levels = sort(filter(<(ef), collect(Iterators.flatten(r.eigenvalues))))
    length(levels) < 2 && return nothing
    jumps = findall(>(gap), diff(levels))
    return isempty(jumps) ? nothing : levels[last(jumps) + 1]
end

"""
    plot_dos(r::DFTKResults; kw...) -> NamedTuple

Open the DOS of `r` in its own window ([`dos_figure`](@ref) with the same
keywords).
"""
function plot_dos(r::DFTKResults; kw...)
    d = dos_figure(r; kw...)
    GLMakie.activate!()
    display(GLMakie.Screen(; title="DOS — " * formula(r.structure)), d.figure)
    return d
end
