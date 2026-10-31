# CIF reader (first data block that has a unit cell).
#
# The asymmetric unit is expanded with the file's symmetry operations
# (`_space_group_symop_operation_xyz` / `_symmetry_equiv_pos_as_xyz`). If the
# file only names the space group, the operations come from spglib's database.

# ---------------------------------------------------------------------------
# Tokenizer and block parser

"""Split CIF text into (token, quoted) pairs."""
function cif_tokens(text::AbstractString)
    tokens = Tuple{String,Bool}[]
    lines = split(text, '\n')
    i = 1
    while i <= length(lines)
        line = rstrip(lines[i], '\r')
        if startswith(line, ';')
            # Multi-line text field, closed by a line starting with ';'.
            buf = [line[2:end]]
            i += 1
            while i <= length(lines) && !startswith(lines[i], ';')
                push!(buf, rstrip(lines[i], '\r'))
                i += 1
            end
            push!(tokens, (strip(join(buf, "\n")), true))
            i += 1
            continue
        end
        pos, last = firstindex(line), lastindex(line)
        while pos <= last
            c = line[pos]
            if isspace(c)
                pos = nextind(line, pos)
            elseif c == '#'
                break
            elseif c == '\'' || c == '"'
                # A quote only closes when followed by whitespace or end of line.
                start = j = nextind(line, pos)
                while j <= last && !(line[j] == c && (j == last || isspace(line[nextind(line, j)])))
                    j = nextind(line, j)
                end
                push!(tokens, (line[start:prevind(line, j)], true))
                pos = nextind(line, j)
            else
                j = pos
                while j <= last && !isspace(line[j])
                    j = nextind(line, j)
                end
                push!(tokens, (line[pos:prevind(line, j)], false))
                pos = j
            end
        end
        i += 1
    end
    return tokens
end

struct CIFBlock
    name::String
    items::Dict{String,String}                  # tag (lowercase) => value
    loops::Vector{Dict{String,Vector{String}}}  # tag (lowercase) => column
end

iskeyword(t, q) = !q && (startswith(t, '_') || lowercase(t) == "loop_" || startswith(lowercase(t), "data_"))

function parse_cif_blocks(text::AbstractString)
    tokens = cif_tokens(text)
    blocks = CIFBlock[]
    k = 1
    n = length(tokens)
    while k <= n
        t, q = tokens[k]
        lt = lowercase(t)
        if !q && startswith(lt, "data_")
            push!(blocks, CIFBlock(t[6:end], Dict(), []))
            k += 1
        elseif isempty(blocks)
            k += 1                                  # anything before the first data_ block
        elseif !q && lt == "loop_"
            k += 1
            tags = String[]
            while k <= n && !tokens[k][2] && startswith(tokens[k][1], '_')
                push!(tags, lowercase(tokens[k][1]))
                k += 1
            end
            values = String[]
            while k <= n && !iskeyword(tokens[k]...)
                push!(values, tokens[k][1])
                k += 1
            end
            isempty(tags) && continue
            length(values) % length(tags) == 0 ||
                error("CIF loop with $(length(tags)) columns has $(length(values)) values")
            cols = Dict(tag => values[j:length(tags):end] for (j, tag) in enumerate(tags))
            push!(blocks[end].loops, cols)
        elseif !q && startswith(t, '_')
            k + 1 <= n || error("CIF tag $t has no value")
            blocks[end].items[lt] = tokens[k+1][1]
            k += 2
        else
            k += 1                                  # stray value; ignore
        end
    end
    return blocks
end

function find_loop(b::CIFBlock, tag)
    for l in b.loops
        haskey(l, tag) && return l
    end
    return nothing
end

"""First of `tags` present as a single item or as a loop column (first value)."""
function cif_item(b::CIFBlock, tags...)
    for t in tags
        haskey(b.items, t) && return b.items[t]
        l = find_loop(b, t)
        l === nothing || return first(l[t])
    end
    return nothing
end

"""Parse a CIF number, dropping a standard uncertainty: `"5.431(2)"` → 5.431. `?` and `.` → nothing."""
function cif_number(s)
    s === nothing && return nothing
    s = strip(s)
    s in ("?", ".") && return nothing
    v = tryparse(Float64, replace(s, r"\(\d+\)$" => ""))
    v === nothing && error("not a number in CIF: '$s'")
    return v
end

# ---------------------------------------------------------------------------
# Symmetry operations

const SymOp = Tuple{SMatrix{3,3,Int,9},Vec3}

parse_fraction(s) = occursin('/', s) ? (/)(parse.(Float64, split(s, '/'))...) : parse(Float64, s)

"""
    parse_symop("-y, x-y, z+1/3") -> (W, w)

Rotation part `W` (integer matrix) and translation `w` of a symmetry
operation in Jones-faithful notation.
"""
function parse_symop(s::AbstractString)
    parts = split(lowercase(replace(s, r"[\s'\"]" => "")), ',')
    length(parts) == 3 || throw(ArgumentError("symmetry operation needs 3 parts: '$s'"))
    W = zeros(Int, 3, 3)
    w = zeros(3)
    for (row, p) in enumerate(parts)
        for m in eachmatch(r"([+-]?)([0-9.]*(?:/[0-9.]+)?)\*?([xyz]?)", p)
            isempty(m.match) && continue
            sgn = m[1] == "-" ? -1 : 1
            num, var = m[2], m[3]
            if isempty(var)
                isempty(num) && continue
                w[row] += sgn * parse_fraction(num)
            else
                coef = isempty(num) ? 1 : round(Int, parse_fraction(num))
                W[row, var[1] - 'x' + 1] += sgn * coef
            end
        end
    end
    return (SMatrix{3,3,Int}(W), Vec3(w))
end

"""Normalize a Hermann–Mauguin symbol for comparison: no spaces, underscores, case."""
normalize_hm(s) = lowercase(replace(s, r"[\s_'\"]" => ""))

"""
    hall_number_from_hm(symbol) -> Int

Find spglib's Hall number for a Hermann–Mauguin symbol such as `"P 6_3 m c"`,
`"F d -3 m :2"` or `"R -3 m H"`. When a group has several settings and none
is given, the first (spglib's default) is used, with a warning.
"""
function hall_number_from_hm(symbol::AbstractString)
    s = strip(symbol)
    choice = ""
    if occursin(':', s)
        s, choice = strip.(split(s, ':'; limit=2))
    else
        # ICSD style: trailing S/Z (origin 1/2) or H/R (axes).
        m = match(r"^(.*\S)\s+([SZHR])$"i, s)
        if m !== nothing
            s = m[1]
            choice = Dict("S" => "1", "Z" => "2", "H" => "H", "R" => "R")[uppercase(m[2])]
        end
    end
    target = normalize_hm(s)
    matches = Int[]
    for h in 1:530
        t = Spglib.get_spacegroup_type(h)
        if target in (normalize_hm(t.international), normalize_hm(t.international_short),
                      normalize_hm(t.international_full))
            push!(matches, h)
        end
    end
    isempty(matches) && error("unknown space group symbol '$symbol'")
    if !isempty(choice)
        for h in matches
            Spglib.get_spacegroup_type(h).choice == choice && return h
        end
        error("space group '$symbol': no setting '$choice'")
    end
    if length(matches) > 1
        choices = [Spglib.get_spacegroup_type(h).choice for h in matches]
        # Rhombohedral groups: CIFs normally use hexagonal axes.
        "H" in choices && return matches[findfirst(==("H"), choices)]
        @warn "Space group '$symbol' has several settings ($(join(choices, ", "))); using '$(choices[1])'. Add e.g. ':2' to the symbol to choose."
    end
    return first(matches)
end

symops_from_hall(h::Integer) =
    [(SMatrix{3,3,Int}(r), Vec3(t)) for (r, t) in zip(Spglib.get_symmetry_from_database(h)...)]

"""Symmetry operations of a CIF block: explicit list, else from the space group name/number."""
function cif_symops(b::CIFBlock, name::AbstractString)
    for tag in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz")
        l = find_loop(b, tag)
        l === nothing || return parse_symop.(l[tag])
        haskey(b.items, tag) && return [parse_symop(b.items[tag])]
    end
    hm = cif_item(b, "_space_group_name_h-m_alt", "_symmetry_space_group_name_h-m")
    if hm !== nothing && !(strip(hm) in ("?", "."))
        return symops_from_hall(hall_number_from_hm(hm))
    end
    num = cif_item(b, "_space_group_it_number", "_symmetry_int_tables_number")
    if num !== nothing && cif_number(num) !== nothing
        n = round(Int, cif_number(num))
        halls = [h for h in 1:530 if Spglib.get_spacegroup_type(h).number == n]
        length(halls) > 1 && @warn "$name: space group $n has several settings; using spglib's default"
        return symops_from_hall(first(halls))
    end
    @warn "$name: no symmetry information; treating as P1"
    return [(SMatrix{3,3,Int}(I), Vec3(0, 0, 0))]
end

# ---------------------------------------------------------------------------
# Reading

"""True if fractional positions f and g are the same site, up to a lattice translation."""
function same_site(L::Lattice, f, g; tol=0.01)
    d = f - g
    return norm(L.matrix * (d - round.(d))) < tol
end

"""Apply all operations to one site and keep the distinct images, wrapped into [0, 1)."""
function orbit(L::Lattice, f0::Vec3, ops::Vector{SymOp})
    images = Vec3[]
    for (W, w) in ops
        f = mod.(W * f0 + w, 1.0)
        f = map(x -> x > 1 - 1e-8 ? 0.0 : x, f)
        any(g -> same_site(L, f, g), images) || push!(images, f)
    end
    return images
end

"""
    read_cif(path) -> Crystal

Reads the first data block that defines a unit cell and expands the
asymmetric unit to the full cell.
"""
read_cif(path::AbstractString) = read_cif(IOBuffer(read(path, String)); name=basename(path))

function read_cif(io::IO; name::AbstractString="CIF")
    blocks = parse_cif_blocks(read(io, String))
    i = findfirst(b -> haskey(b.items, "_cell_length_a"), blocks)
    i === nothing && error("$name: no data block with a unit cell (_cell_length_a)")
    b = blocks[i]

    cell = [cif_number(get(b.items, t, nothing)) for t in
            ("_cell_length_a", "_cell_length_b", "_cell_length_c",
             "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")]
    any(isnothing, cell[1:3]) && error("$name: incomplete cell lengths")
    cell[4:6] = [something(x, 90.0) for x in cell[4:6]]
    L = Lattice(cell...)

    sites = find_loop(b, "_atom_site_fract_x")
    sites === nothing && error("$name: no _atom_site_fract_x loop (Cartesian-only CIFs are not supported)")
    labels = get(sites, "_atom_site_label", nothing)
    types = get(sites, "_atom_site_type_symbol", labels)
    types === nothing && error("$name: atom sites have neither _atom_site_type_symbol nor _atom_site_label")
    occ = get(sites, "_atom_site_occupancy", nothing)

    ops = cif_symops(b, name)
    species = Symbol[]
    frac = Vec3[]
    partial = String[]
    for k in eachindex(sites["_atom_site_fract_x"])
        f0 = Vec3(cif_number(sites[c][k]) for c in ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"))
        sym = normalize_symbol(types[k])
        o = occ === nothing ? 1.0 : something(cif_number(occ[k]), 1.0)
        o < 0.999 && push!(partial, labels === nothing ? string(sym) : labels[k])
        for f in orbit(L, f0, ops)
            push!(species, sym)
            push!(frac, f)
        end
    end
    isempty(partial) ||
        @warn "$name: partially occupied sites $(join(partial, ", ")) are drawn as fully occupied"
    if any(same_site(L, frac[j], frac[k]) for j in eachindex(frac) for k in 1:j-1)
        @warn "$name: two sites share a position (disorder?); both are kept"
    end

    title = something(cif_item(b, "_chemical_name_mineral", "_chemical_name_systematic",
                                "_chemical_formula_sum"), b.name)
    return Crystal(L, species, frac; title=strip(title))
end
