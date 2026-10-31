# Unit-cell box: 12 edges. The three edges leaving the origin are colored
# like the a/b/c axes and labelled, as in Materials Studio.

const CELL_EDGES = let corners = [(i, j, k) for i in 0:1 for j in 0:1 for k in 0:1]
    # Pairs of corners that differ in exactly one coordinate.
    [(p, q) for p in corners for q in corners if p < q && sum(abs.(p .- q)) == 1]
end

"""Draw the cell box; returns the plots created."""
function draw_cell!(ax, L::Lattice, vs::ViewSettings)
    corner(p) = Point3f(frac_to_cart(L, Vec3(p)))
    origin_edges = [((0, 0, 0), (1, 0, 0)), ((0, 0, 0), (0, 1, 0)), ((0, 0, 0), (0, 0, 1))]

    others = Point3f[]
    for (p, q) in CELL_EDGES
        (p, q) in origin_edges && continue
        push!(others, corner(p), corner(q))
    end
    plots = Any[Makie.linesegments!(ax, others; color=vs.cell_color, linewidth=vs.cell_linewidth)]

    axes_pts = Point3f[]
    for (p, q) in origin_edges
        push!(axes_pts, corner(p), corner(q))
    end
    axcolors = [vs.axis_colors[i] for i in 1:3 for _ in 1:2]
    push!(plots, Makie.linesegments!(ax, axes_pts; color=axcolors, linewidth=2vs.cell_linewidth))

    if vs.show_axis_labels
        a, b, c = cell_vectors(L)
        tips = [Point3f(v + 0.6normalize(v)) for v in (a, b, c)]   # 0.6 Å past the edge
        push!(plots, Makie.text!(ax, tips; text=["a", "b", "c"], color=collect(vs.axis_colors),
            fontsize=20, font=:bold, align=(:center, :center)))
    end
    return plots
end

"""
Draw the box of a molecule (a structure that is not periodic) as thin dashed
grey lines: it is the cell written to input files, but only a container.
"""
function draw_molecule_box!(ax, L::Lattice, vs::ViewSettings)
    corner(p) = Point3f(frac_to_cart(L, Vec3(p)))
    pts = Point3f[]
    for (p, q) in CELL_EDGES
        push!(pts, corner(p), corner(q))
    end
    return Any[Makie.linesegments!(ax, pts; color=vs.molecule_box_color, linewidth=vs.cell_linewidth,
                                   linestyle=:dash)]
end
