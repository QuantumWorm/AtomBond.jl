# VASP POSCAR / CONTCAR reader and writer (VASP 5 format with a species line;
# VASP 4 files are accepted if the species can be taken from the comment line).

"""
    read_poscar(path) -> Crystal
"""
read_poscar(path::AbstractString) = open(io -> read_poscar(io; name=basename(path)), path)

function read_poscar(io::IO; name::AbstractString="POSCAR")
    lines = readlines(io)   # positional format: keep blank lines (the comment may be empty)
    length(lines) >= 7 || error("$name: too short to be a POSCAR file")
    title = strip(lines[1])

    scale = parse(Float64, split(lines[2])[1])
    vecs = [Vec3(parse.(Float64, split(lines[i])[1:3])) for i in 3:5]
    M = hcat(vecs...)
    # A negative scale factor means "target volume" in VASP.
    factor = scale < 0 ? cbrt(-scale / abs(det(M))) : scale
    M = M * factor

    k = 6
    tokens = split(lines[k])
    if all(t -> tryparse(Int, t) !== nothing, tokens)
        # VASP 4: no species line; try the comment line.
        counts = parse.(Int, tokens)
        names = split(title)
        length(names) >= length(counts) ||
            error("$name: no species line and the comment line does not list $(length(counts)) element names")
        species_names = names[1:length(counts)]
    else
        species_names = tokens
        k += 1
        counts = parse.(Int, split(lines[k]))
        length(counts) == length(species_names) ||
            error("$name: $(length(species_names)) species but $(length(counts)) counts")
    end
    k += 1

    selective = lowercase(strip(lines[k]))[1] == 's'   # Selective dynamics: T/F flags per coordinate
    selective && (k += 1)
    mode = lowercase(strip(lines[k]))[1]
    cartesian = mode in ('c', 'k')
    k += 1

    n = sum(counts)
    length(lines) >= k + n - 1 || error("$name: expected $n coordinate lines")
    species = reduce(vcat, [fill(normalize_symbol(s), c) for (s, c) in zip(species_names, counts)])

    L = Lattice(M)
    frac = map(0:n-1) do i
        v = Vec3(parse.(Float64, split(lines[k+i])[1:3]))
        cartesian ? cart_to_frac(L, v * factor) : v
    end
    # The flags always refer to fractional coordinates; F = fixed.
    fixed = !selective ? nothing : map(0:n-1) do i
        w = split(lines[k+i])
        length(w) >= 6 ? Tuple(uppercase(t[1]) == 'F' for t in w[4:6]) : FREE_ATOM
    end
    Crystal(L, species, frac; title=title, fixed)
end

"""
    write_poscar(path, crystal; cartesian=false)

Write a VASP 5 POSCAR. Atoms are grouped by species (in first-appearance order).
"""
write_poscar(path::AbstractString, c::Crystal; kwargs...) = open(io -> write_poscar(io, c; kwargs...), path, "w")

function write_poscar(io::IO, c::Crystal; cartesian::Bool=false)
    println(io, isempty(c.title) ? formula(c) : c.title)
    println(io, "1.0")
    for v in cell_vectors(c.lattice)
        @printf(io, "  %20.14f %20.14f %20.14f\n", v...)
    end
    spc = species_list(c)
    println(io, "  ", join(string.(spc), "  "))
    println(io, "  ", join([count(==(s), c.species) for s in spc], "  "))
    # Fixed atoms: Selective dynamics, F = fixed (flags refer to fractional coordinates).
    selective = any(any, c.fixed)
    selective && println(io, "Selective dynamics")
    println(io, cartesian ? "Cartesian" : "Direct")
    for s in spc, i in eachindex(c.species)
        c.species[i] == s || continue
        v = cartesian ? frac_to_cart(c.lattice, c.frac[i]) : c.frac[i]
        @printf(io, "  %20.14f %20.14f %20.14f", v...)
        selective && print(io, "   ", join((x ? "F" : "T" for x in c.fixed[i]), " "))
        println(io)
    end
end
