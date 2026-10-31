# Editing operations. Each returns a new Crystal and leaves the input alone,
# so the viewer's undo history is simply a stack of crystals.

function check_indices(c::Crystal, idx)
    isempty(idx) && throw(ArgumentError("no atoms given"))
    bad = filter(i -> !(1 <= i <= natoms(c)), collect(idx))
    isempty(bad) || throw(ArgumentError("no atom(s) $(join(bad, ", ")) (structure has $(natoms(c)))"))
    return unique(idx)
end

"""
    known_element(s) -> Symbol

Parse an element symbol typed by the user (`"ge"`, `"Ge"`); throws for
anything that is not in the element table.
"""
function known_element(s)
    str = strip(string(s))
    sym = isempty(str) ? :X : Symbol(uppercase(str[1:1]) * lowercase(str[2:end]))
    haskey(ELEMENTS, sym) || throw(ArgumentError("unknown element '$str'"))
    return sym
end

"""
    delete_atoms(c, indices) -> Crystal
"""
function delete_atoms(c::Crystal, idx)
    idx = check_indices(c, idx)
    keep = setdiff(1:natoms(c), idx)
    Crystal(c.lattice, c.species[keep], c.frac[keep]; pbc=c.pbc, title=c.title, fixed=c.fixed[keep], charges=c.charges[keep],
            magmoms=c.magmoms[keep],
            cell_constraints=c.cell_constraints)
end

"""
    substitute(c, indices, element) -> Crystal

Change the element of the given atoms (e.g. a Ge substitution in Si).
"""
function substitute(c::Crystal, idx, el)
    idx = check_indices(c, idx)
    sym = known_element(el)
    species = copy(c.species)
    species[idx] .= sym
    # A new element on a site: its old charge and moment no longer apply.
    charges = [s == c.species[i] ? c.charges[i] : NaN for (i, s) in enumerate(species)]
    magmoms = [s == c.species[i] ? c.magmoms[i] : NaN for (i, s) in enumerate(species)]
    Crystal(c.lattice, species, c.frac; pbc=c.pbc, title=c.title, fixed=c.fixed, cell_constraints=c.cell_constraints,
            charges, magmoms)
end

"""
    add_atom(c, element, position; cartesian=false, min_distance=0.5) -> Crystal

Add one atom at a fractional (or, with `cartesian=true`, Cartesian Å)
position. Refuses positions closer than `min_distance` Å to an existing
atom, including periodic images.
"""
function add_atom(c::Crystal, el, pos; cartesian::Bool=false, min_distance::Real=0.5)
    sym = known_element(el)
    f = cartesian ? cart_to_frac(c.lattice, pos) : Vec3(pos)
    trial = Crystal(c.lattice, [c.species; sym], [c.frac; [f]]; pbc=c.pbc, title=c.title,
                    fixed=[c.fixed; [FREE_ATOM]], cell_constraints=c.cell_constraints, charges=[c.charges; NaN],
                    magmoms=[c.magmoms; NaN])
    n = natoms(trial)
    for j in 1:n-1
        d = min_image_distance(trial, n, j)
        d < min_distance && throw(ArgumentError(@sprintf("too close to %s (%.2f Å)", atom_label(c, j), d)))
    end
    return trial
end

"""
    move_atoms(c, indices, shift) -> Crystal

Translate the given atoms by a Cartesian vector (Å).
"""
function move_atoms(c::Crystal, idx, shift)
    idx = check_indices(c, idx)
    df = c.lattice.matrix \ Vec3(shift)
    frac = copy(c.frac)
    frac[idx] .+= Ref(df)
    with_positions(c, frac)
end

"""
    translate_atoms(c, vector; fractional=false, atoms=all) -> Crystal

Shift all atoms (or the indices in `atoms`) by `vector`: Cartesian Å, or
fractional coordinates of a, b, c with `fractional=true`. Atoms may end up
outside the cell; [`wrap`](@ref) maps them back.
"""
function translate_atoms(c::Crystal, vector; fractional::Bool=false, atoms=nothing)
    length(vector) == 3 || throw(ArgumentError("a shift needs three numbers"))
    df = fractional ? Vec3(vector) : cart_to_frac(c.lattice, Vec3(vector))
    idx = atoms === nothing ? (1:natoms(c)) : check_indices(c, atoms)
    frac = copy(c.frac)
    frac[idx] .+= Ref(df)
    with_positions(c, frac)
end

"""
    rotation_matrix(axis, angle) -> 3×3 matrix

Rotation by `angle` degrees about the Cartesian `axis` (right-hand rule:
counterclockwise when the axis points toward the viewer).
"""
function rotation_matrix(axis, angle::Real)
    n = norm(axis)
    n > 0 || throw(ArgumentError("the rotation axis is zero"))
    u = Vec3(axis) / n
    K = SMatrix{3,3,Float64}(0, u[3], -u[2], -u[3], 0, u[1], u[2], -u[1], 0)
    θ = deg2rad(angle)
    return one(SMatrix{3,3,Float64}) + sin(θ) * K + (1 - cos(θ)) * K * K
end

"""
    rotate_atoms(c, angle, axis; center=centroid, atoms=all, images=nothing) -> Crystal

Rotate all atoms (or the indices in `atoms`) by `angle` degrees about the
Cartesian `axis` through `center` (Å; default: the centroid of the rotated
atoms). The cell stays fixed. `images` (one lattice shift per rotated atom)
rotates those periodic copies instead, so a group that crosses a cell face
turns as one piece.
"""
function rotate_atoms(c::Crystal, angle::Real, axis; center=nothing, atoms=nothing, images=nothing)
    idx = atoms === nothing ? collect(1:natoms(c)) : collect(atoms)
    isempty(idx) && throw(ArgumentError("no atoms given"))
    allunique(idx) || throw(ArgumentError("atoms are listed twice"))
    check_indices(c, idx)
    img = images === nothing ? fill(Shift(0, 0, 0), length(idx)) : Shift.(images)
    length(img) == length(idx) || throw(ArgumentError("need one image per rotated atom"))
    r = [frac_to_cart(c.lattice, c.frac[i] + s) for (i, s) in zip(idx, img)]
    p = center === nothing ? sum(r) / length(r) : Vec3(center)
    R = rotation_matrix(axis, angle)
    frac = copy(c.frac)
    for (k, i) in enumerate(idx)
        frac[i] = cart_to_frac(c.lattice, p + R * (r[k] - p))
    end
    with_positions(c, frac)
end

"""
    closest_contact(c, atoms; limit=200_000) -> (i, j, distance) or nothing

The closest pair between `atoms` and any other atom (periodic images
included), e.g. to spot overlaps after moving a group. `nothing` when there
is no other atom, or when checking would take more than `limit` pairs.
"""
function closest_contact(c::Crystal, atoms; limit::Integer=200_000)
    n = natoms(c)
    (n < 2 || length(atoms) * n > limit) && return nothing
    best = nothing
    for i in atoms, j in 1:n
        j == i && continue
        d = min_image_distance(c, i, j)
        (best === nothing || d < best[3]) && (best = (i, j, d))
    end
    return best
end

"""
    set_position(c, index, position; cartesian=false, min_distance=0.5) -> Crystal

Put atom `index` at a fractional (or, with `cartesian=true`, Cartesian Å)
position. Refuses positions closer than `min_distance` Å to another atom,
including periodic images.
"""
function set_position(c::Crystal, i::Integer, pos; cartesian::Bool=false, min_distance::Real=0.5)
    i = only(check_indices(c, [i]))
    f = cartesian ? cart_to_frac(c.lattice, Vec3(pos)) : Vec3(pos)
    frac = copy(c.frac)
    frac[i] = f
    trial = with_positions(c, frac)
    for j in 1:natoms(c)
        j == i && continue
        d = min_image_distance(trial, i, j)
        d < min_distance && throw(ArgumentError(@sprintf("too close to %s (%.2f Å)", atom_label(c, j), d)))
    end
    return trial
end

"""
    replace_species(c, from, to) -> Crystal

Change every atom of element `from` into `to` (e.g. all Ga → In).
"""
function replace_species(c::Crystal, from, to)
    src, dst = known_element(from), known_element(to)
    idx = findall(==(src), c.species)
    isempty(idx) && throw(ArgumentError("there is no $src in $(formula(c))"))
    src == dst && throw(ArgumentError("$src is already $dst"))
    return substitute(c, idx, dst)
end
