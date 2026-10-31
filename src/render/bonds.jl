# Bonds as instanced cylinders: one meshscatter call for all half-bonds.
# The unit cylinder runs from z = 0 to z = 1 with radius 1; each instance is
# rotated onto the bond direction and scaled to (radius, radius, length).
# The :line style draws the same half-bonds as colored line segments.

halfbond_color(dm::DisplayModel, h::HalfBond, vs::ViewSettings) =
    vs.bond_color === nothing ? element(dm.crystal.species[h.index]).color : vs.bond_color

function draw_bonds!(ax, dm::DisplayModel, vs::ViewSettings)
    isempty(dm.halfbonds) && return nothing
    vs.style === :line && return draw_bond_lines!(ax, dm, vs)

    starts = Point3f[]
    dirs = Vec3f[]
    scales = Vec3f[]
    colors = RGBf[]
    r = Float32(vs.bond_radius)
    for h in dm.halfbonds
        v = h.to - h.from
        len = norm(v)
        len > 0 || continue
        push!(starts, Point3f(h.from))
        push!(dirs, Vec3f(v / len))
        push!(scales, Vec3f(r, r, len))
        push!(colors, halfbond_color(dm, h, vs))
    end
    cylinder = Makie.Tessellation(Cylinder(Point3f(0, 0, 0), Point3f(0, 0, 1), 1.0f0), vs.bond_resolution)
    Makie.meshscatter!(ax, starts;
        marker=cylinder, rotation=dirs, markersize=scales, color=colors,
        specular=vs.specular, shininess=vs.shininess, ssao=vs.ssao)
end

function draw_bond_lines!(ax, dm::DisplayModel, vs::ViewSettings)
    pts = Point3f[]
    colors = RGBf[]
    for h in dm.halfbonds
        c = halfbond_color(dm, h, vs)
        push!(pts, Point3f(h.from), Point3f(h.to))
        push!(colors, c, c)
    end
    Makie.linesegments!(ax, pts; color=colors, linewidth=vs.line_width)
end
