# Camera: framing, projection and "view along" presets.
#
# The LScene gets one Camera3D (`setup_camera!`); everything afterwards only
# changes its settings, so redrawing plots never moves the view.

"""Create the 3D camera once. `center=false` stops Makie re-centering when plots are added."""
setup_camera!(lscene, vs::ViewSettings) =
    Makie.cam3d!(lscene.scene; center=false, fov=vs.fov,
                 projectiontype=vs.orthographic ? Makie.Orthographic : Makie.Perspective)

camera(lscene) = Makie.cameracontrols(lscene.scene)

"""
Center and radius of a sphere enclosing the cell and all drawn atoms and
bonds (for a molecule: the atoms, plus its box when the cell is shown).
"""
function bounding_sphere(dm::DisplayModel, vs::ViewSettings)
    c = dm.crystal
    pts = [a.position for a in dm.atoms]
    append!(pts, (h.to for h in dm.halfbonds))   # stubs reach past the outer atoms
    if any(c.pbc) || vs.show_cell || isempty(pts)          # a shown box (also a molecule's) is framed too
        append!(pts, vec([frac_to_cart(c.lattice, Vec3(p...)) for p in Iterators.product(0:1, 0:1, 0:1)]))
        center = frac_to_cart(c.lattice, Vec3(0.5, 0.5, 0.5))
    else
        center = sum(pts) / length(pts)
    end
    pad = maximum((atom_radius(s, vs) for s in c.species); init=0.0)
    radius = maximum(norm(p - center) for p in pts) + pad
    return center, radius
end

"""
    look_along!(lscene, dm, settings, direction; up)

Place the camera on the `direction` side of the structure (direction points
from the structure toward the eye), far enough away that everything fits.
"""
function look_along!(lscene, dm::DisplayModel, vs::ViewSettings, direction::AbstractVector;
                     up::AbstractVector=Vec3(0, 0, 1))
    center, radius = bounding_sphere(dm, vs)
    dir = normalize(Vec3(direction))
    up = Vec3(up)
    # An up vector parallel to the view direction is undefined; pick another.
    if norm(cross(dir, normalize(up))) < 1e-3
        up = abs(dir[3]) < 0.9 ? Vec3(0, 0, 1) : Vec3(0, 1, 0)
    end
    cam = camera(lscene)
    set_projection!(cam, vs)
    # Distance at which the bounding sphere just fills the view. Makie's
    # orthographic camera shows ±distance vertically, independent of fov.
    dist = vs.orthographic ? 1.05radius : radius / sind(vs.fov / 2)
    Makie.update_cam!(lscene.scene, cam, Vec3f(center + dist * dir), Vec3f(center), Vec3f(up))
    return cam
end

"""
Switch between perspective and orthographic projection. Makie's orthographic
zoom moves the eye toward the structure, so its near clip plane would cut
into the atoms; fixed planes far on either side of the eye avoid that.
"""
function set_projection!(cam, vs::ViewSettings)
    cam.fov[] = vs.fov
    if vs.orthographic
        cam.settings.projectiontype[] = Makie.Orthographic
        cam.settings.clipping_mode[] = :static
        cam.near[] = -1000.0   # Å, relative to the eye
        cam.far[] = 1000.0
    else
        cam.settings.projectiontype[] = Makie.Perspective
        cam.settings.clipping_mode[] = :adaptive
        cam.near[] = 0.1       # × eye distance
        cam.far[] = 100.0
    end
end

"""Re-frame along the current viewing direction (e.g. after a projection change)."""
function reframe!(lscene, dm::DisplayModel, vs::ViewSettings)
    cam = camera(lscene)
    look_along!(lscene, dm, vs, cam.eyeposition[] - cam.lookat[]; up=cam.upvector[])
end

"""Default oblique view (`settings.view_direction`)."""
reset_camera!(lscene, dm::DisplayModel, vs::ViewSettings) =
    look_along!(lscene, dm, vs, vs.view_direction)

"""
    view_direction(lattice, which) -> (direction, up)

Presets: `:a, :b, :c` look down a lattice vector; `:astar, :bstar, :cstar`
look down a reciprocal vector, i.e. straight onto the (100), (010), (001)
planes. A 3-vector is taken as a lattice direction [uvw].
"""
function view_direction(L::Lattice, which)
    a, b, c = cell_vectors(L)
    recip = inv(L.matrix)'          # columns ∝ a*, b*, c*
    which === :a && return a, c
    which === :b && return b, c
    which === :c && return c, b
    which === :astar && return recip[:, 1], c
    which === :bstar && return recip[:, 2], c
    which === :cstar && return recip[:, 3], b
    which isa AbstractVector && return L.matrix * Vec3(which), c
    throw(ArgumentError("unknown view preset $which"))
end

"""
    rotate_view!(lscene, angle)

Turn the picture by `angle` degrees in the screen plane (counter-clockwise
for a positive angle) about the line of sight: the viewing direction stays
the same, only the up direction of the camera changes. Useful after a view
along an axis, to put another axis horizontal or vertical.
"""
function rotate_view!(lscene, angle::Real)
    cam = camera(lscene)
    eye, lookat = Vec3(cam.eyeposition[]), Vec3(cam.lookat[])
    d = normalize(eye - lookat)                         # toward the viewer
    up = Vec3(cam.upvector[])
    up = up - dot(up, d) * d                            # the part in the screen plane
    norm(up) < 1e-9 && (up = abs(d[3]) < 0.9 ? Vec3(0, 0, 1) - d[3] * d : Vec3(0, 1, 0) - d[2] * d)
    up = normalize(up)
    # The picture turns counter-clockwise when the camera's up vector turns
    # clockwise about the line of sight (Rodrigues rotation by -angle about d).
    θ = -deg2rad(angle)
    new = up * cos(θ) + cross(d, up) * sin(θ)
    Makie.update_cam!(lscene.scene, cam, Vec3f(eye), Vec3f(lookat), Vec3f(new))
    return cam
end

"""Apply a preset from [`view_direction`](@ref)."""
function view_along!(lscene, dm::DisplayModel, vs::ViewSettings, which)
    dir, up = view_direction(dm.crystal.lattice, which)
    look_along!(lscene, dm, vs, dir; up)
end
