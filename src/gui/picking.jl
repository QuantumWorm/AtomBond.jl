# Selecting atoms with the mouse and measuring between them.
#
#   click            select an atom (replaces the selection)
#   Shift + click    add an atom to the selection, or remove it again
#   click on empty   clear the selection          (Esc does the same)
#   right-click      menu for the selected atoms (contextmenu.jl)
#
# Selecting 2, 3 or 4 atoms shows their distance, angle or dihedral, in the
# order they were clicked. A left-drag still rotates the view: only a press
# and release without movement counts as a click.

const SELECTION_HINT = "Click an atom to select it; Shift+click adds atoms (2 = distance, 3 = angle, 4 = dihedral)."
const CLICK_SLOP = 4.0   # pixels the mouse may move and still count as a click

"""
Drawn atoms for the current selection keys. A key whose image is no longer
drawn (the atom was moved or wrapped) moves to the atom's first drawn copy;
keys of atoms that are gone are dropped.
"""
function selected_atoms(v::Viewer)
    lookup = Dict((a.index, a.image) => a for a in v.display.atoms)
    first_copy = Dict{Int,Shift}()
    foreach(a -> get!(first_copy, a.index, a.image), v.display.atoms)
    map!(k -> haskey(lookup, k) || !haskey(first_copy, k[1]) ? k : (k[1], first_copy[k[1]]), v.selection, v.selection)
    unique!(v.selection)
    filter!(k -> haskey(lookup, k), v.selection)
    return [lookup[k] for k in v.selection]
end

"""Redraw the highlight/measurement overlay and the info line."""
function update_selection!(v::Viewer)
    close_context_menu!(v)      # its title and items describe the old selection
    foreach(p -> delete!(v.lscene.scene, p), v.overlay)
    atoms = selected_atoms(v)
    v.overlay = draw_selection!(v.lscene, v.display, v.settings, atoms)
    v.info[] = selection_info(v.display, atoms)
    return v
end

"""Name an atom; images outside the home cell get a prime: Si1′."""
display_label(c::Crystal, a::DisplayAtom) =
    atom_label(c, a.index) * (iszero(a.image) ? "" : "′")

"""" fixed x y" for an atom's fixed fractional coordinates, or ""."""
fixed_text(f) = any(f) ? "   fixed " * (all(f) ? "x y z" : join(("x", "y", "z")[collect(f)], " ")) : ""

"""One line describing the selection: atom details, or the measurement."""
function selection_info(dm::DisplayModel, atoms::Vector{DisplayAtom})
    c = dm.crystal
    n = length(atoms)
    n == 0 && return SELECTION_HINT
    if n == 1
        a = atoms[1]
        f = c.frac[a.index] + a.image
        nb = count(b -> b.i == a.index || b.j == a.index, dm.bonds)
        return @sprintf("%s   frac (%.4f, %.4f, %.4f)   cart (%.4f, %.4f, %.4f) Å   %d bond%s",
                        display_label(c, a), f..., a.position..., nb, nb == 1 ? "" : "s") *
               fixed_text(c.fixed[a.index]) *
               (isnan(c.charges[a.index]) ? "" : "   charge $(format_charge(c.charges[a.index])) e") *
               (isnan(c.magmoms[a.index]) ? "" : "   moment $(format_magmom(c.magmoms[a.index]))")
    end
    names = join((display_label(c, a) for a in atoms), " – ")
    m = measure([a.position for a in atoms])
    m === nothing && return "$n atoms selected: $names"
    kind, value = m
    return "$names   $(kind) $(format_measure(kind, value))"
end

"""
    select!(viewer, index, image=(0,0,0); add=false)

Select an atom programmatically (as a click would). `add=true` acts like
Shift+click: toggles the atom in the current selection.
"""
function select!(v::Viewer, index::Integer, image=Shift(0, 0, 0); add::Bool=false)
    key = (Int(index), Shift(image))
    if add
        k = findfirst(==(key), v.selection)
        k === nothing ? push!(v.selection, key) : deleteat!(v.selection, k)
    else
        v.selection = [key]
    end
    update_selection!(v)
end

function clear_selection!(v::Viewer)
    empty!(v.selection)
    update_selection!(v)
end

# Pixel offsets tried around the cursor: the measurement lines are drawn on
# top and pass through atom centres, so the exact pixel may be a line.
const PICK_OFFSETS = [(0, 0); [(round(Int, r * cos(t)), round(Int, r * sin(t)))
                               for r in (3, 6) for t in range(0, 2π; length=9)[1:8]]]

"""The atom at window pixel `xy` (via the atom plot or a selection halo), or `nothing`."""
function atom_at(v::Viewer, xy)
    plot, idx = Makie.pick(v.lscene.scene, xy)
    plot === nothing && return nothing
    if plot === v.atom_plot
        atoms = drawn_atoms(v.display, v.settings)
    elseif !isempty(v.overlay) && plot === first(v.overlay)
        atoms = selected_atoms(v)      # the halo around a selected atom
    else
        return nothing                 # bonds, cell, labels, measurement
    end
    return 1 <= idx <= length(atoms) ? atoms[idx] : nothing
end

"""What is under the mouse: a DisplayAtom, or `nothing`."""
function atom_under_mouse(v::Viewer)
    xy = Vec2f(Makie.events(v.lscene.scene).mouseposition[])
    for off in PICK_OFFSETS
        a = atom_at(v, xy .+ Vec2f(off))
        a === nothing || return a
    end
    return nothing
end

function install_picking!(v::Viewer)
    scene = v.lscene.scene
    ev = Makie.events(scene)
    press = Dict{Makie.Mouse.Button,Any}()
    # Priority above the camera's handlers, but never consume the event, so
    # rotation (left drag) and panning (right drag) keep working.
    Makie.on(ev.mousebutton; priority=10) do e
        e.button in (Makie.Mouse.left, Makie.Mouse.right) || return Makie.Consume(false)
        if e.action == Makie.Mouse.press
            press[e.button] = Makie.is_mouseinside(scene) ? ev.mouseposition[] : nothing
        elseif e.action == Makie.Mouse.release && get(press, e.button, nothing) !== nothing
            moved = norm(Vec2f(ev.mouseposition[]) - Vec2f(press[e.button]))
            press[e.button] = nothing
            moved <= CLICK_SLOP && (e.button == Makie.Mouse.left ? handle_click!(v) : handle_right_click!(v))
        end
        return Makie.Consume(false)
    end
    Makie.on(ev.keyboardbutton; priority=10) do e
        if e.key == Makie.Keyboard.escape && e.action == Makie.Keyboard.press
            clear_selection!(v)
        end
        return Makie.Consume(false)
    end
end

function handle_click!(v::Viewer)
    a = atom_under_mouse(v)
    click_ripple!(v, Makie.events(v.figure.scene).mouseposition[], a === nothing ? :empty : :atom)
    a === nothing || click_sound()
    shift = Makie.ispressed(v.lscene.scene, Makie.Keyboard.left_shift | Makie.Keyboard.right_shift)
    if a === nothing
        shift || clear_selection!(v)
    else
        select!(v, a.index, a.image; add=shift)
    end
end
