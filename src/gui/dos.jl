# Density-of-states plot of a finished DFTK run (its own window).

const DOS_FILL = RGBAf(0.20, 0.45, 0.75, 0.30)    # filled states
const DOS_LINE = RGBf(0.12, 0.30, 0.55)
const DOS_DOWN_FILL = RGBAf(0.80, 0.35, 0.20, 0.30)
const DOS_DOWN_LINE = RGBf(0.60, 0.22, 0.10)

"""Ways to split the DOS into projected curves (menu label, `by` of [`pdos`](@ref))."""
const PDOS_GROUPINGS = [("No projection", :none), ("Element and orbital", :element_orbital), ("Element", :element),
                        ("Orbital (s, p, d)", :orbital), ("Each atom", :atom)]

"""The projection choices for `r` (each atom only for up to 12 atoms)."""
pdos_choices(r::DFTResults) = r.projections === nothing ? PDOS_GROUPINGS[1:1] :
                               filter(o -> o[2] !== :atom || natoms(r.structure) <= 12, PDOS_GROUPINGS)

"""
Draw the projected DOS of `r` grouped `by` on `ax` (energies relative to
`ef`; `horizontal`: DOS along x), replacing the earlier curves in `plots`
and the legend in `legend`.
"""
function draw_pdos!(ax, plots, legend, r::DFTResults, by::Symbol, width::Real, ef; horizontal::Bool=false)
    foreach(p -> delete!(ax, p), plots)
    empty!(plots)
    legend[] === nothing || delete!(legend[])
    legend[] = nothing
    if by !== :none
        pd = pdos(r; by, width)
        colors = Makie.wong_colors()
        for (i, (n, m)) in enumerate(pd.curves), σ in axes(m, 2)
            x, y = pd.energies .- ef, (σ == 1 ? 1 : -1) .* m[:, σ]
            push!(plots, Makie.lines!(ax, horizontal ? y : x, horizontal ? x : y; color=colors[mod1(i, length(colors))],
                                      linewidth=1.4, label=σ == 1 ? n : nothing))
        end
    end
    if by !== :none || (r.n_spin == 2 && !horizontal)
        legend[] = Makie.axislegend(ax; position=:rt, framevisible=false, labelsize=12, rowgap=0,
                                    backgroundcolor=RGBAf(1, 1, 1, 0.85))
    end
    return plots
end

"""
    dos_figure(r::DFTResults; width=0.1, title="", dir=nothing, name="dos") -> (; figure, axis, width, save, ...)

A figure with the density of states of `r` ([`dos`](@ref)) against the
energy relative to the Fermi level: the filled states (below the Fermi
level) shaded, the projected DOS ([`pdos`](@ref)) as coloured curves when the
run has orbital projections (*Projected DOS* menu; `projection` is the
initial grouping, `:element_orbital` by default), the Fermi level dashed, spin down mirrored below the axis
for a spin-polarized run. A slider sets the Gaussian broadening (eV) and the
plot follows it. Deep levels (semicore states, separated from the valence
levels by an empty range of more than `deep_gap` eV, 5 by default) are left
out of the view unless *Show deep levels* is on; the CSV always has them. With `dir`, a Save button writes `name.dos.png`
and `name.dos.csv` there. The figure is not displayed; see [`plot_dos`](@ref).
"""
function dos_figure(r::DFTResults; width::Real=0.1, title::AbstractString="", dir=nothing, name::AbstractString="dos",
                    deep_gap::Real=5, projection::Symbol=:element_orbital)
    isempty(r.eigenvalues) && throw(ArgumentError("no levels in these results (a relaxation keeps none)"))
    fig = Makie.Figure(; size=(760, 520), backgroundcolor=:white)
    ef = something(r.fermi, 0.0)
    head = !isempty(title) ? title : !isempty(r.structure.title) ? r.structure.title : formula(r.structure)
    nk = length(first(dos_levels(r))) ÷ r.n_spin
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
    choices = pdos_choices(r)
    proj = Makie.Observable(:none)
    if length(choices) > 1
        gp = Makie.GridLayout(fig[row, 1:3]; tellwidth=false, halign=:left)
        Makie.Label(gp[1, 1], "Projected DOS"; halign=:left)
        start = something(findfirst(o -> o[2] === projection, choices), 1)
        menu = Makie.Menu(gp[1, 2]; options=choices, default=choices[start][1], width=200)
        Makie.connect!(proj, menu.selection)
        proj[] = choices[start][2]
        row += 1
    end
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
    pplots, legend = Any[], Ref{Any}(nothing)
    redraw_pdos(_...) = draw_pdos!(ax, pplots, legend, r, proj[], slider.value[], ef)
    Makie.on(redraw_pdos, proj)
    Makie.on(redraw_pdos, slider.value)
    redraw_pdos()
    function fit!(_...)                   # limits from the data in view (autolimits needs a displayed figure)
        d = data[]
        xs = d.energies .- ef
        lo = showdeep[] || valence_bottom === nothing ? first(xs) : valence_bottom - ef - 4 * slider.value[]
        hi = min(last(xs), 10.0)                     # the highest computed bands are not reliable
        top = maximum(d.dos[(xs .>= lo) .& (xs .<= hi), :]) * 1.08
        Makie.limits!(ax, lo, hi, r.n_spin == 2 ? -top : 0.0, top)
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
    return (; figure=fig, axis=ax, width=slider, save, data, showdeep, projection=proj, pdos_plots=pplots)
end

"""
Lowest valence level (eV) when deep (semicore) levels lie below it, separated
by an empty range wider than `gap` eV below the Fermi level; otherwise `nothing`.
"""
function deep_levels_edge(r::DFTResults; gap::Real=5)
    ef = something(r.fermi, Inf)
    levels = sort(filter(<(ef), collect(Iterators.flatten(r.eigenvalues))))
    length(levels) < 2 && return nothing
    jumps = findall(>(gap), diff(levels))
    return isempty(jumps) ? nothing : levels[last(jumps) + 1]
end

"""
    plot_dos(r::DFTResults; kw...) -> NamedTuple

Open the DOS of `r` in its own window ([`dos_figure`](@ref) with the same
keywords).
"""
function plot_dos(r::DFTResults; kw...)
    d = dos_figure(r; kw...)
    GLMakie.activate!()
    display(GLMakie.Screen(; title="DOS — " * formula(r.structure)), d.figure)
    return d
end

"""
    write_bands(path, r) -> path

The band structure of `r` ([`BandPath`](@ref)) as CSV: distance along the
path (1/Å), the label of the point (if any), then every band (eV, relative to
the Fermi level), spin up first for a spin-polarized run.
"""
function write_bands(path::AbstractString, r::DFTResults)
    p = r.path
    p === nothing && throw(ArgumentError("no band structure in these results"))
    ef = something(r.fermi, 0.0)
    cols = [r.n_spin == 1 ? "band$b" : "$(σ == 1 ? "up" : "down")$b" for (σ, m) in enumerate(p.eigenvalues) for b in axes(m, 2)]
    open(path, "w") do io
        println(io, "k_distance_inv_A,label,", join(cols, ","))
        for k in eachindex(p.distances)
            i = findfirst(t -> isapprox(t, p.distances[k]; atol=1e-6), p.ticks)
            label = i === nothing ? "" : replace(p.labels[i], "," => " ")
            vals = (@sprintf("%.5f", m[k, b] - ef) for m in p.eigenvalues for b in axes(m, 2))
            println(io, @sprintf("%.6f", p.distances[k]), ",", label, ",", join(vals, ","))
        end
    end
    return path
end

"""Polyline through every band of `m` (k × band), NaN between bands and path pieces, energies minus `ef`."""
function band_lines(p::BandPath, m::AbstractMatrix, ef)
    pieces = [s:(i < length(p.branch_starts) ? p.branch_starts[i+1] - 1 : length(p.distances))
              for (i, s) in enumerate(p.branch_starts)]
    xs, ys = Float64[], Float64[]
    for b in axes(m, 2), piece in pieces
        append!(xs, p.distances[piece]); push!(xs, NaN)
        append!(ys, m[piece, b] .- ef); push!(ys, NaN)
    end
    return xs, ys
end

"""
    bands_figure(r::DFTResults; width=0.1, window=nothing, title="", dir=nothing, name="bands") -> NamedTuple

The band structure of `r` along the high-symmetry path (left) next to the
DOS (right, sideways) on a shared energy axis relative to the Fermi level
(dashed); `window` is the energy range shown (eV, relative to the Fermi
level), by default from the bottom of the valence levels (deep semicore
levels left out) to 8 eV. The valence-band maximum
and conduction-band minimum are marked and the gap is given in the subtitle
([`band_edges`](@ref)). A slider sets the DOS broadening; with `dir`, a Save
button writes `name.bands.png`, `name.bands.csv` and `name.dos.csv` there.
The DOS panel shows the projected DOS ([`pdos`](@ref)) as coloured curves when
the run has orbital projections (*Projected DOS* menu, `projection`).
Not displayed; see [`plot_bands`](@ref).
"""
function bands_figure(r::DFTResults; width::Real=0.1, window=nothing, title::AbstractString="", dir=nothing,
                      name::AbstractString="bands", projection::Symbol=:element_orbital)
    p = r.path
    p === nothing && throw(ArgumentError("no band structure in these results (run with bands=true)"))
    ef = something(r.fermi, 0.0)
    head = !isempty(title) ? title : !isempty(r.structure.title) ? r.structure.title : formula(r.structure)
    e = band_edges(r)
    gapnote = is_metallic(r) ? "metallic" : e === nothing ? "" :
              e.direct ? @sprintf("direct gap %.3f eV at %s", e.gap, e.vbm_at) :
              @sprintf("indirect gap %.3f eV (VBM %s, CBM %s)", e.gap, e.vbm_at, e.cbm_at)
    fig = Makie.Figure(; size=(900, 560), backgroundcolor=:white)
    heading = Makie.GridLayout(fig[0, 1:2])
    Makie.Label(heading[1, 1], head; font=:bold, fontsize=16)
    Makie.Label(heading[2, 1], "DFTK $(r.version)" * (isempty(gapnote) ? "" : ", " * gapnote) *
                               (r.n_spin == 2 ? ", spin up blue, down red" : ""); fontsize=13, color=:gray30)
    Makie.rowgap!(heading, 2)
    ylabel = r.fermi === nothing ? "E (eV)" : "E − E_F (eV)"
    axb = Makie.Axis(fig[1, 1]; ylabel, xticks=(p.ticks, p.labels), xgridvisible=false, ygridvisible=false)
    axd = Makie.Axis(fig[1, 2]; xlabel="DOS (states/eV)", yticklabelsvisible=false, xgridvisible=false,
                     ygridvisible=false)
    Makie.colsize!(fig.layout, 2, Makie.Relative(0.25))
    Makie.linkyaxes!(axb, axd)
    Makie.vlines!(axb, p.ticks; color=:gray70, linewidth=0.8)
    for (σ, m) in enumerate(p.eigenvalues)
        xs, ys = band_lines(p, m, ef)
        Makie.lines!(axb, xs, ys; color=σ == 1 ? DOS_LINE : DOS_DOWN_LINE, linewidth=1.3)
    end
    r.fermi === nothing || Makie.hlines!(axb, [0.0]; color=:gray30, linestyle=:dash, linewidth=1)
    if e !== nothing && e.direct !== nothing
        kv = argmin(k -> abs(maximum(m[k, b] for m in p.eigenvalues for b in axes(m, 2) if m[k, b] <= ef; init=-Inf) - e.vbm),
                    eachindex(p.distances))
        kc = argmin(k -> abs(minimum(m[k, b] for m in p.eigenvalues for b in axes(m, 2) if m[k, b] > ef; init=Inf) - e.cbm),
                    eachindex(p.distances))
        Makie.scatter!(axb, [p.distances[kv], p.distances[kc]], [e.vbm - ef, e.cbm - ef];
                       color=[:black, :black], markersize=8)
    end
    slider = Makie.Slider(fig[2, 1]; range=0.02:0.01:0.5, startvalue=clamp(width, 0.02, 0.5))
    wlabel = Makie.lift(w -> @sprintf("DOS broadening %.2f eV", w), slider.value)
    Makie.Label(fig[2, 2], wlabel; halign=:left)
    data = Makie.lift(w -> dos(r; width=w), slider.value)
    y = Makie.lift(d -> d.energies .- ef, data)
    for σ in 1:r.n_spin
        sgn = σ == 1 ? 1 : -1
        x = Makie.lift(d -> sgn .* d.dos[:, σ], data)
        filled = Makie.lift((ys, xs) -> [yi <= 0 ? xi : 0.0 for (yi, xi) in zip(ys, xs)], y, x)
        Makie.band!(axd, y, Makie.lift(v -> zero(v), filled), filled; direction=:y,
                    color=σ == 1 ? DOS_FILL : DOS_DOWN_FILL)
        Makie.lines!(axd, x, y; color=σ == 1 ? DOS_LINE : DOS_DOWN_LINE, linewidth=1.2)
    end
    r.fermi === nothing || Makie.hlines!(axd, [0.0]; color=:gray30, linestyle=:dash, linewidth=1)
    allb = reduce(vcat, [vec(m) for m in p.eigenvalues]) .- ef
    if window === nothing
        deep = deep_levels_edge(r)
        bottom = deep === nothing ? minimum(allb) : deep - ef
        lo, hi = bottom - 1.0, min(8.0, maximum(allb))
    else
        lo, hi = Float64.(window)
    end
    Makie.xlims!(axb, first(p.distances), last(p.distances))
    function fit!(d)
        sel = (d.energies .- ef .>= lo) .& (d.energies .- ef .<= hi)
        top = maximum(d.dos[sel, :]; init=1.0) * 1.08
        Makie.limits!(axd, r.n_spin == 2 ? -top : 0.0, top, lo, hi)
        Makie.ylims!(axb, lo, hi)
    end
    Makie.on(fit!, data)
    fit!(data[])
    choices = pdos_choices(r)
    proj = Makie.Observable(:none)
    row = 3
    if length(choices) > 1
        gp = Makie.GridLayout(fig[row, 1:2]; tellwidth=false, halign=:left)
        Makie.Label(gp[1, 1], "Projected DOS"; halign=:left)
        start = something(findfirst(o -> o[2] === projection, choices), 1)
        menu = Makie.Menu(gp[1, 2]; options=choices, default=choices[start][1], width=200)
        Makie.connect!(proj, menu.selection)
        proj[] = choices[start][2]
        row += 1
    end
    pplots, legend = Any[], Ref{Any}(nothing)
    redraw_pdos(_...) = draw_pdos!(axd, pplots, legend, r, proj[], slider.value[], ef; horizontal=true)
    Makie.on(redraw_pdos, proj)
    Makie.on(redraw_pdos, slider.value)
    redraw_pdos()
    save = dir === nothing ? nothing : Makie.Button(fig[row, 1:2]; label="Save PNG and CSV to " * basename(dir) * "/")
    if save !== nothing
        Makie.on(save.clicks) do _
            write_bands(joinpath(dir, name * ".bands.csv"), r)
            write_dos(joinpath(dir, name * ".dos.csv"), r; width=slider.value[])
            Makie.save(joinpath(dir, name * ".bands.png"), fig; update=false)
            save.label[] = "Saved $name.bands.png, $name.bands.csv and $name.dos.csv"
        end
    end
    return (; figure=fig, bands_axis=axb, dos_axis=axd, width=slider, save, data, edges=e, projection=proj,
            pdos_plots=pplots)
end

"""
    plot_bands(r::DFTResults; kw...) -> NamedTuple

Open the band structure and DOS of `r` in their own window
([`bands_figure`](@ref) with the same keywords).
"""
function plot_bands(r::DFTResults; kw...)
    d = bands_figure(r; kw...)
    GLMakie.activate!()
    display(GLMakie.Screen(; title="Bands — " * formula(r.structure)), d.figure)
    return d
end
