# VASP output: XDATCAR (MD / relaxation trajectory) and vasprun.xml.

"""
    read_xdatcar(path or io; frame=:last) -> Crystal

One frame of an XDATCAR (by default the last). Variable-cell runs, which
repeat the header before every frame, are handled.
"""
read_xdatcar(path::AbstractString; frame=:last) =
    open(io -> read_xdatcar(io; frame, title=basename(path)), path)

function read_xdatcar(io::IO; frame=:last, title::AbstractString="XDATCAR")
    lines = readlines(io)
    configs = findall(l -> occursin(r"^\s*(Direct|Cartesian)\s+configuration"i, l), lines)
    isempty(configs) && throw(ArgumentError("no 'Direct configuration=' lines"))
    k = frame === :last ? configs[end] : configs[frame]
    # The header for this frame: the nearest one before it (only the first,
    # for a fixed cell). A header is title, scale, 3 vectors, species, counts.
    isheader(j) = j >= 7 && tryparse(Float64, strip(lines[j-5])) !== nothing &&
                  all(t -> tryparse(Int, t) !== nothing, split(lines[j]))
    h = findlast(j -> isheader(j), 1:k-1)
    h === nothing && throw(ArgumentError("no header before frame"))
    scale = parse(Float64, split(lines[h-5])[1])
    M = hcat((Vec3(parse.(Float64, split(lines[h-5+i])[1:3])) for i in 1:3)...) * scale
    names, counts = split(lines[h-1]), parse.(Int, split(lines[h]))
    species = reduce(vcat, [fill(normalize_symbol(s), n) for (s, n) in zip(names, counts)])
    L = Lattice(M)
    xyz = [Vec3(parse.(Float64, split(lines[k+i])[1:3])) for i in 1:sum(counts)]
    step = something(match(r"=\s*(\d+)", lines[k]), (; captures=["?"])).captures[1]
    t = "$(strip(lines[h-6])), step $step"
    occursin(r"^\s*Cartesian"i, lines[k]) && return crystal_from_cartesian(L, species, [scale * r for r in xyz]; title=t)
    return Crystal(L, species, xyz; title=t)
end

"""Rows of numbers of `<varray name="name">` inside `xml`."""
function vasprun_varray(xml::AbstractString, name::AbstractString)
    m = match(Regex("<varray\\s+name=\"$name\"\\s*>(.*?)</varray>", "s"), xml)
    m === nothing && return nothing
    return [Vec3(parse.(Float64, split(v[1]))) for v in eachmatch(r"<v[^>]*>([^<]*)</v>", m[1])]
end

"""
    read_vasprun(path or io) -> Crystal

The final structure of a vasprun.xml (`finalpos`, or the last structure
written when the run did not finish).
"""
read_vasprun(path::AbstractString) = open(io -> read_vasprun(io; title=basename(path)), path)

function read_vasprun(io::IO; title::AbstractString="vasprun.xml")
    xml = read(io, String)
    atoms = match(r"<array\s+name=\"atoms\"\s*>(.*?)</array>"s, xml)
    atoms === nothing && throw(ArgumentError("no atoms array (not a vasprun.xml?)"))
    species = [normalize_symbol(m[1]) for m in eachmatch(r"<rc>\s*<c>\s*([^<]+?)\s*</c>", atoms[1])]
    final = match(r"<structure\s+name=\"finalpos\"\s*>(.*?)</structure>"s, xml)
    s = if final !== nothing
        final[1]
    else
        all_s = collect(eachmatch(r"<structure[^>]*>(.*?)</structure>"s, xml))
        isempty(all_s) && throw(ArgumentError("no structure"))
        all_s[end][1]
    end
    basis, pos = vasprun_varray(s, "basis"), vasprun_varray(s, "positions")
    (basis === nothing || pos === nothing) && throw(ArgumentError("structure without basis or positions"))
    length(pos) == length(species) || throw(ArgumentError("$(length(species)) atoms but $(length(pos)) positions"))
    return Crystal(Lattice(basis...), species, pos; title=final === nothing ? "$title (last step, unfinished)" : title)
end
