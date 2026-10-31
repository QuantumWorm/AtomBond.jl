# Small a/b/c axes indicator in the lower-left corner of the 3D view. It is a
# separate scene with its own camera that copies the main camera's orientation,
# so it turns with the structure but never moves or zooms.
#
# It hangs off the figure's root scene, not the LScene: GLMakie applies SSAO
# only to scenes without children, so a child would switch SSAO off.

const TRIAD_SIZE = 110   # pixels

"""
Add the triad. `lattice` is an Observable so the arrows follow cell changes
(e.g. switching to the primitive cell).
"""
function add_triad!(root, lscene, lattice::Makie.Observable{Lattice}, vs::ViewSettings)
    viewport = Makie.lift(vp -> Makie.Rect2i(Makie.origin(vp) .+ 6, (TRIAD_SIZE, TRIAD_SIZE)),
                          lscene.scene.viewport)
    triad = Makie.Scene(root; viewport, clear=false)
    # No mouse/keyboard control of its own: it only follows the main camera.
    tcam = Makie.cam3d!(triad; center=false, fov=30,
        mouse_rotationspeed=0, mouse_translationspeed=0, mouse_zoomspeed=0,
        keyboard_rotationspeed=0, keyboard_translationspeed=0, keyboard_zoomspeed=0)

    dirs = Makie.lift(L -> [Vec3f(normalize(v)) for v in cell_vectors(L)], lattice)
    colors = collect(vs.axis_colors)
    Makie.arrows3d!(triad, fill(Point3f(0), 3), dirs; color=colors,
        shaftradius=0.05, tipradius=0.12, tiplength=0.3, markerscale=1)
    Makie.text!(triad, Makie.lift(ds -> [Point3f(1.35d) for d in ds], dirs); text=["a", "b", "c"],
        color=colors, fontsize=16, font=:bold, align=(:center, :center))

    main = camera(lscene)
    Makie.onany(main.eyeposition, main.lookat, main.upvector) do eye, lookat, up
        dir = normalize(Vec3f(eye - lookat))
        Makie.update_cam!(triad, tcam, 6.5f0 * dir, Vec3f(0), Vec3f(up))
    end
    Makie.notify(main.eyeposition)   # sync once now
    return triad
end
