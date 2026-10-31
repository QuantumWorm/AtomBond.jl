# What to draw for a crystal, worked out without touching any graphics.
#
# - Atoms are drawn where they are, also outside the cell (as in Materials
#   Studio; Build → Wrap atoms into cell moves them back). With
#   `wrap_atoms=true` they are drawn wrapped into the cell instead.
# - Along periodic directions an atom sitting on a cell face (frac ≈ 0 or 1)
#   is also drawn on the opposite face, so a corner atom appears at all 8
#   corners, as in Materials Studio.
# - Bonds between two drawn atoms are drawn in full, split at the midpoint
#   and colored by each end's element.
# - A bond to an atom that is not drawn (it lies outside the cell) becomes a
#   half-bond stub up to the midpoint, if `bond_stubs` is on.

"""A drawn atom: crystal index `index`, displaced by `image` lattice vectors."""
struct DisplayAtom
    index::Int
    image::Shift
    position::Vec3
end

"""One colored cylinder: half of a bond, from `from` to `to`, belonging to atom `index`."""
struct HalfBond
    index::Int
    from::Vec3
    to::Vec3
end

struct DisplayModel
    crystal::Crystal            # wrapped copy actually drawn
    atoms::Vector{DisplayAtom}
    halfbonds::Vector{HalfBond}
    bonds::Vector{Bond}
end

"""
Wrap into [0, 1) along periodic directions, snapping values within `tol` of 1
down to 0 so they are treated as boundary atoms.
"""
function wrap_for_display(c::Crystal, tol::Real)
    frac = map(c.frac) do f
        Vec3(ntuple(3) do k
            c.pbc[k] || return f[k]
            x = mod(f[k], 1.0)
            x > 1 - tol ? x - 1 : x
        end)
    end
    with_positions(c, frac)
end

"""
Whole cells between atom `f0` and its wrapped copy `f`, i.e. where to draw it.
Along a direction where the atom is on a face it is drawn on both faces of
the cell it is next to (both faces of the unit cell for 0 or 1).
"""
function home_image(f0, f, pbc, tol)
    Shift(ntuple(3) do k
        pbc[k] || return 0
        m = round(Int, f0[k] - f[k])
        abs(f[k]) < tol && m >= 1 ? m - 1 : m
    end)
end

function boundary_images(c::Crystal, tol::Real, home)
    atoms = DisplayAtom[]
    for (i, f) in enumerate(c.frac)
        # Along each periodic direction where the atom is on the face, it also
        # appears one cell over.
        options = ntuple(k -> (c.pbc[k] && abs(f[k]) < tol) ? (0, 1) : (0,), 3)
        for s1 in options[1], s2 in options[2], s3 in options[3]
            s = Shift(s1, s2, s3) + home[i]
            push!(atoms, DisplayAtom(i, s, frac_to_cart(c.lattice, f + s)))
        end
    end
    return atoms
end

"""
    build_display(crystal, settings) -> DisplayModel
"""
function build_display(c0::Crystal, vs::ViewSettings)
    c = wrap_for_display(c0, vs.boundary_tol)
    # Bonds are found in the wrapped copy; images count from it.
    home = vs.wrap_atoms ? fill(Shift(0, 0, 0), natoms(c)) :
           [home_image(f0, f, c.pbc, vs.boundary_tol) for (f0, f) in zip(c0.frac, c.frac)]
    atoms = vs.show_boundary_atoms ? boundary_images(c, vs.boundary_tol, home) :
            [DisplayAtom(i, home[i], frac_to_cart(c.lattice, f + home[i])) for (i, f) in enumerate(c.frac)]

    bonds = bonds_visible(vs) ? find_bonds(c; tolerance=vs.bond_tolerance, cutoffs=vs.bond_cutoffs) : Bond[]
    halfbonds = HalfBond[]
    isempty(bonds) && return DisplayModel(c, atoms, halfbonds, bonds)

    # Neighbors of each atom as (j, shift), from both ends of every bond.
    adj = [Tuple{Int,Shift}[] for _ in 1:natoms(c)]
    for b in bonds
        push!(adj[b.i], (b.j, b.shift))
        push!(adj[b.j], (b.i, -b.shift))
    end

    drawn = Set((a.index, a.image) for a in atoms)
    pos(i, s) = frac_to_cart(c.lattice, c.frac[i] + s)
    # Every drawn atom draws its own half of each of its bonds. A bond between
    # two drawn atoms therefore comes out whole and two-colored; a bond to an
    # undrawn atom is just the stub.
    for a in atoms, (j, s) in adj[a.index]
        other = (j, a.image + s)
        (other in drawn || vs.bond_stubs) || continue
        mid = (a.position + pos(other...)) / 2
        push!(halfbonds, HalfBond(a.index, a.position, mid))
    end
    return DisplayModel(c, atoms, halfbonds, bonds)
end
