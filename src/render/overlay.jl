# Things drawn on top of the structure: atom labels, the selection highlight
# and the measurement (distance / angle / dihedral) of the selected atoms.

# Labels and measurements are drawn with `overdraw` so atoms never hide them;
# a white outline keeps the text readable over bonds.
const TEXT_STYLE = (overdraw=true, strokecolor=:white, strokewidth=2)

function draw_labels!(ax, dm::DisplayModel, vs::ViewSettings)
    atoms = drawn_atoms(dm, vs)
    isempty(atoms) && return nothing
    Makie.text!(ax, [Point3f(a.position) for a in atoms];
        text=[atom_label(dm.crystal, a.index) for a in atoms],
        fontsize=vs.label_size, color=:black, align=(:left, :bottom), offset=(4, 4),
        TEXT_STYLE...)
end

"""
Markers on fixed atoms: × when all three fractional coordinates are fixed,
+ when only some are. Drawn on top so they stay visible.
"""
function draw_fixed!(ax, dm::DisplayModel, vs::ViewSettings)
    c = dm.crystal
    atoms = [a for a in drawn_atoms(dm, vs) if any(c.fixed[a.index])]
    isempty(atoms) && return nothing
    Makie.scatter!(ax, [Point3f(a.position) for a in atoms];
        marker=[all(c.fixed[a.index]) ? :xcross : :cross for a in atoms], markersize=13,
        color=:black, strokecolor=:white, strokewidth=1.5, overdraw=true)
end

"""Bonded neighbours of each atom: (j, shift, length)."""
function bond_neighbours(dm::DisplayModel)
    adj = [Tuple{Int,Shift,Float64}[] for _ in 1:natoms(dm.crystal)]
    for b in dm.bonds
        push!(adj[b.i], (b.j, b.shift, b.length))
        push!(adj[b.j], (b.i, -b.shift, b.length))
    end
    return adj
end

"""
Positions and values of bond angles: for each drawn atom, the angle between
every pair of its bonds, placed inside the angle (along the bisector, 40 %
of the shorter bond out). `between` limits it to angles whose three atoms
are all in that set of (index, image) keys; otherwise every bond counts,
including stubs when those are shown.
"""
function bond_angle_labels(dm::DisplayModel, vs::ViewSettings; between=nothing)
    c = dm.crystal
    adj = bond_neighbours(dm)
    drawn = Set((a.index, a.image) for a in dm.atoms)
    pts, texts = Point3f[], String[]
    for a in dm.atoms
        between === nothing || (a.index, a.image) in between || continue
        nb = Vec3[]
        for (j, s, _) in adj[a.index]
            other = (j, a.image + s)
            between === nothing ? (other in drawn || vs.bond_stubs || continue) : (other in between || continue)
            push!(nb, frac_to_cart(c.lattice, c.frac[j] + other[2]) - a.position)
        end
        for x in 1:length(nb), y in x+1:length(nb)
            u, w = nb[x], nb[y]
            bis = normalize(u) + normalize(w)
            # A straight angle has no bisector: put the label beside the atom.
            dir = norm(bis) > 1e-3 ? normalize(bis) : normalize(cross(u, abs(u[3]) < 0.9norm(u) ? Vec3(0, 0, 1) : Vec3(1, 0, 0)))
            push!(pts, Point3f(a.position + 0.4 * min(norm(u), norm(w)) * dir))
            push!(texts, @sprintf("%.1f°", bond_angle(a.position + u, a.position, a.position + w)))
        end
    end
    return pts, texts
end

"""Bond-angle labels (°); see [`bond_angle_labels`](@ref)."""
function draw_bond_angles!(ax, dm::DisplayModel, vs::ViewSettings; between=nothing)
    pts, texts = bond_angle_labels(dm, vs; between)
    isempty(pts) && return nothing
    Makie.text!(ax, pts; text=texts, fontsize=vs.label_size - 1, font=:bold, color=vs.bond_angle_color,
                align=(:center, :center), TEXT_STYLE...)
end

"""
Positions and values of dihedral angles i–j–k–l: for each bond j–k (once),
every bonded neighbour i of j and l of k, placed around the bond (see
[`dihedral_label_position`](@ref)). Torsions through a straight angle (undefined) are left out. `between`
limits it to dihedrals whose four atoms are all in that set of keys;
otherwise the central bond must join two drawn atoms, and the outer atoms
may be outside the cell when stubs are shown.
"""
function dihedral_labels(dm::DisplayModel, vs::ViewSettings; between=nothing)
    c = dm.crystal
    adj = bond_neighbours(dm)
    drawn = Set((a.index, a.image) for a in dm.atoms)
    pos(key) = frac_to_cart(c.lattice, c.frac[key[1]] + key[2])
    usable(key) = between === nothing ? (key in drawn || vs.bond_stubs) : key in between
    neighbours(key) = [(j, key[2] + s) for (j, s, _) in adj[key[1]]]
    pts, texts, seen = Point3f[], String[], Set{NTuple{3,Int}}()
    for a in dm.atoms
        jkey = (a.index, a.image)
        between === nothing || jkey in between || continue
        for kkey in neighbours(jkey)
            (between === nothing ? kkey in drawn : kkey in between) || continue
            mid = Tuple(round.(Int, 1000 .* (pos(jkey) + pos(kkey)) ./ 2))    # each bond once
            mid in seen && continue
            push!(seen, mid)
            pj, pk = pos(jkey), pos(kkey)
            for ikey in neighbours(jkey), lkey in neighbours(kkey)
                (ikey == kkey || lkey == jkey || ikey == lkey) && continue
                (usable(ikey) && usable(lkey)) || continue
                pi, pl = pos(ikey), pos(lkey)
                (bond_angle(pi, pj, pk) > 179 || bond_angle(pj, pk, pl) > 179) && continue
                push!(pts, Point3f(dihedral_label_position(pi, pj, pk, pl)))
                φ = dihedral(pi, pj, pk, pl)
                push!(texts, @sprintf("%.1f°", φ < -179.95 ? φ + 360 : φ))      # -180 and 180 alike
            end
        end
    end
    return pts, texts
end

"""
Where to put the label of dihedral i–j–k–l: on a ring around the middle of
bond j–k, towards the outer atoms i and l (their directions seen along the
bond, averaged), 40 % of the bond length out. For an anti (180°) pair those
directions cancel, so the label goes at right angles to atom i, 25 % out,
which keeps the several dihedrals of one bond apart.
"""
function dihedral_label_position(pi, pj, pk, pl)
    L = norm(pk - pj)
    axis = (pk - pj) / L
    across(v) = (w = v - dot(v, axis) * axis; normalize(w))
    ui, ul = across(pi - pj), across(pl - pk)
    b = ui + ul
    mid = (pj + pk) / 2
    norm(b) > 0.2 && return mid + 0.4L * normalize(b)
    return mid + 0.25L * normalize(cross(axis, ui))
end

"""Dihedral-angle labels (°); see [`dihedral_labels`](@ref)."""
function draw_dihedrals!(ax, dm::DisplayModel, vs::ViewSettings; between=nothing)
    pts, texts = dihedral_labels(dm, vs; between)
    isempty(pts) && return nothing
    Makie.text!(ax, pts; text=texts, fontsize=vs.label_size - 1, font=:bold, color=vs.dihedral_color,
                align=(:center, :center), TEXT_STYLE...)
end

"""Coordinate system used for coordinate labels: `:fractional` or `:cartesian`."""
coordinate_system(c::Crystal, vs::ViewSettings) =
    vs.coordinates === :auto ? (any(c.pbc) ? :fractional : :cartesian) : vs.coordinates

"""
Positions and texts of coordinate labels for drawn atoms: "(x, y, z)" with 3
decimals, fractional or Cartesian Å (see [`coordinate_system`](@ref)). A
copy of an atom on a cell face shows its own coordinates (e.g. 1.000).
"""
function coordinate_labels(dm::DisplayModel, vs::ViewSettings, atoms)
    c = dm.crystal
    frac = coordinate_system(c, vs) === :fractional
    texts = map(atoms) do a
        x = frac ? cart_to_frac(c.lattice, a.position) : a.position
        @sprintf("(%.3f, %.3f, %.3f)", (abs(t) < 5e-4 ? 0.0 : t for t in x)...)
    end
    return [Point3f(a.position) for a in atoms], texts
end

"""Coordinate labels below-right of the atoms (atom names go above-right)."""
function draw_coordinates!(ax, dm::DisplayModel, vs::ViewSettings, atoms)
    isempty(atoms) && return nothing
    pts, texts = coordinate_labels(dm, vs, atoms)
    Makie.text!(ax, pts; text=texts, fontsize=vs.label_size - 1, color=vs.coordinate_color,
                align=(:left, :top), offset=(4, -4), TEXT_STYLE...)
end

"""Charge text: "+0.46", "−0.46" (true minus sign), "0.00"."""
format_charge(q) = (s = @sprintf("%+.2f", q); s == "+0.00" || s == "-0.00" ? "0.00" : replace(s, "-" => "−"))

"""
Positions, texts and colours of charge labels for the drawn atoms whose
charge is known.
"""
function charge_labels(dm::DisplayModel, vs::ViewSettings, atoms)
    c = dm.crystal
    known = [a for a in atoms if !isnan(c.charges[a.index])]
    q = [c.charges[a.index] for a in known]
    colors = [x < 0 ? vs.negative_charge_color : vs.positive_charge_color for x in q]
    return [Point3f(a.position) for a in known], format_charge.(q), colors
end

"""Charge labels above-left of the atoms (names go above-right, coordinates below-right)."""
function draw_charges!(ax, dm::DisplayModel, vs::ViewSettings, atoms)
    pts, texts, colors = charge_labels(dm, vs, atoms)
    isempty(pts) && return nothing
    Makie.text!(ax, pts; text=texts, color=colors, fontsize=vs.label_size, font=:bold,
                align=(:right, :bottom), offset=(-4, 4), TEXT_STYLE...)
end

"""Moment text: "+2.15 μB", "−2.15 μB", "0.00 μB"."""
format_magmom(m) = format_charge(m) * " μB"

"""
Positions, texts and colours of magnetic-moment labels for the drawn atoms
whose moment is known (positive orange, negative teal).
"""
function magmom_labels(dm::DisplayModel, vs::ViewSettings, atoms)
    c = dm.crystal
    known = [a for a in atoms if !isnan(c.magmoms[a.index])]
    m = [c.magmoms[a.index] for a in known]
    colors = [x < 0 ? vs.spin_down_color : vs.spin_up_color for x in m]
    return [Point3f(a.position) for a in known], format_magmom.(m), colors
end

"""Moment labels below-left of the atoms (charges go above-left, coordinates below-right)."""
function draw_magmoms!(ax, dm::DisplayModel, vs::ViewSettings, atoms)
    pts, texts, colors = magmom_labels(dm, vs, atoms)
    isempty(pts) && return nothing
    Makie.text!(ax, pts; text=texts, color=colors, fontsize=vs.label_size, font=:bold,
                align=(:right, :top), offset=(-4, -4), TEXT_STYLE...)
end

"""
Midpoints and lengths of the drawn bonds, each bond once. With `between`
(a set of (index, image) keys of drawn atoms) only bonds with both ends in
it; otherwise every bond of a drawn atom, including stubs to atoms outside
the cell when those are shown.
"""
function bond_length_labels(dm::DisplayModel, vs::ViewSettings; between=nothing)
    c = dm.crystal
    adj = bond_neighbours(dm)
    drawn = Set((a.index, a.image) for a in dm.atoms)
    pts, texts, seen = Point3f[], String[], Set{NTuple{3,Int}}()
    for a in dm.atoms, (j, s, d) in adj[a.index]
        other = (j, a.image + s)
        if between === nothing
            other in drawn || vs.bond_stubs || continue
        else
            ((a.index, a.image) in between && other in between) || continue
        end
        mid = (a.position + frac_to_cart(c.lattice, c.frac[j] + other[2])) / 2
        key = Tuple(round.(Int, 1000 .* mid))           # each bond once (both ends find it)
        key in seen && continue
        push!(seen, key)
        push!(pts, Point3f(mid))
        push!(texts, @sprintf("%.3f", d))
    end
    return pts, texts
end

"""Bond-length labels (Å) at bond midpoints; see [`bond_length_labels`](@ref)."""
function draw_bond_lengths!(ax, dm::DisplayModel, vs::ViewSettings; between=nothing)
    pts, texts = bond_length_labels(dm, vs; between)
    isempty(pts) && return nothing
    Makie.text!(ax, pts; text=texts, fontsize=vs.label_size, font=:bold, color=vs.bond_length_color,
                align=(:center, :center), TEXT_STYLE...)
end

"""Text for a measurement value."""
format_measure(kind, value) =
    kind === :distance ? @sprintf("%.4f Å", value) : @sprintf("%.2f°", value)

"""
    draw_selection!(ax, dm, settings, selected) -> plots

Halo around each selected atom; for 2–4 atoms also dashed lines joining them
in selection order and the measured value.
"""
function draw_selection!(ax, dm::DisplayModel, vs::ViewSettings, selected::Vector{DisplayAtom})
    plots = Any[]
    isempty(selected) && return plots
    c = dm.crystal
    pts = [Point3f(a.position) for a in selected]
    radii = Float32[max(1.35 * atom_radius(c.species[a.index], vs), 0.35) for a in selected]
    # Plain alpha blending (not GLMakie's order-independent transparency, which
    # is composited over everything and would fade the text drawn on top).
    push!(plots, Makie.meshscatter!(ax, pts; markersize=radii, color=vs.highlight_color,
                                    transparency=false, shading=Makie.NoShading))
    # Coordinates of the selected atoms,
    if vs.bond_lengths === :selected && vs.show_coordinates
        p = draw_coordinates!(ax, dm, vs, selected)
        p === nothing || push!(plots, p)
    end
    # charges,
    if vs.bond_lengths === :selected && vs.show_charges
        p = draw_charges!(ax, dm, vs, selected)
        p === nothing || push!(plots, p)
    end
    # magnetic moments,
    if vs.bond_lengths === :selected && vs.show_magmoms
        p = draw_magmoms!(ax, dm, vs, selected)
        p === nothing || push!(plots, p)
    end
    # lengths of the bonds among 3 or more selected atoms (2 show the measured distance),
    # and of the bond angles among 4 or more (3 show the measured angle).
    if vs.bond_lengths === :selected && bonds_visible(vs) && length(selected) >= 3
        keys = Set((a.index, a.image) for a in selected)
        if vs.show_bond_lengths
            p = draw_bond_lengths!(ax, dm, vs; between=keys)
            p === nothing || push!(plots, p)
        end
        if vs.show_bond_angles && length(selected) >= 4
            p = draw_bond_angles!(ax, dm, vs; between=keys)
            p === nothing || push!(plots, p)
        end
        if vs.show_dihedrals && length(selected) >= 5             # 4 show the measured dihedral
            p = draw_dihedrals!(ax, dm, vs; between=keys)
            p === nothing || push!(plots, p)
        end
    end

    m = measure([a.position for a in selected])
    m === nothing && return plots
    kind, value = m
    segs = Point3f[]
    for k in 1:length(pts)-1
        push!(segs, pts[k], pts[k+1])
    end
    push!(plots, Makie.linesegments!(ax, segs; color=vs.measure_color, linewidth=2,
                                     linestyle=:dash, overdraw=true))
    # Distance: midpoint. Angle: the vertex. Dihedral: middle of the central bond.
    at = kind === :distance ? (pts[1] + pts[2]) / 2 :
         kind === :angle ? pts[2] : (pts[2] + pts[3]) / 2
    push!(plots, Makie.text!(ax, [at]; text=[format_measure(kind, value)], fontsize=16,
                             font=:bold, color=vs.measure_color, align=(:left, :bottom),
                             offset=(14, 10), TEXT_STYLE...))
    return plots
end
