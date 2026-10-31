# Quantum ESPRESSO: pw.x input files (any common ibrav) and pw.x output
# (the last structure of a relax / vc-relax, or the input structure of an
# scf run). The pw.x writer is in export.jl.

"""Namelist variables of a pw.x input: Dict(lowercase name => value text)."""
function qe_namelists(text::AbstractString)
    vars = Dict{String,String}()
    for m in eachmatch(r"&\w+(.*?)\n\s*/"s, text)
        body = replace(m.captures[1], r"!.*" => "")
        for a in eachmatch(r"([A-Za-z_][\w()]*)\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s]+)", body)
            vars[lowercase(a.captures[1])] = strip(a.captures[2], ['\'', '"'])
        end
    end
    return vars
end

qe_float(s) = parse(Float64, replace(lowercase(s), 'd' => 'e'))

"""Cell vectors (Å) of Bravais lattice `ibrav` from celldm (celldm[1] = a in Bohr)."""
function qe_lattice(ibrav::Integer, celldm::Vector{Float64})
    a = celldm[1] * BOHR
    ba, ca = celldm[2], celldm[3]
    b, c = a * ba, a * ca
    V = Vec3
    vs = if ibrav == 1
        (V(a, 0, 0), V(0, a, 0), V(0, 0, a))
    elseif ibrav == 2
        (a / 2 * V(-1, 0, 1), a / 2 * V(0, 1, 1), a / 2 * V(-1, 1, 0))
    elseif ibrav == 3
        (a / 2 * V(1, 1, 1), a / 2 * V(-1, 1, 1), a / 2 * V(-1, -1, 1))
    elseif ibrav == -3
        (a / 2 * V(-1, 1, 1), a / 2 * V(1, -1, 1), a / 2 * V(1, 1, -1))
    elseif ibrav == 4
        (V(a, 0, 0), V(-a / 2, a * sqrt(3) / 2, 0), V(0, 0, c))
    elseif ibrav == 5
        cα = celldm[4]
        tx, ty, tz = sqrt((1 - cα) / 2), sqrt((1 - cα) / 6), sqrt((1 + 2cα) / 3)
        (a * V(tx, -ty, tz), a * V(0, 2ty, tz), a * V(-tx, -ty, tz))
    elseif ibrav == 6
        (V(a, 0, 0), V(0, a, 0), V(0, 0, c))
    elseif ibrav == 7
        (a / 2 * V(1, -1, ca), a / 2 * V(1, 1, ca), a / 2 * V(-1, -1, ca))
    elseif ibrav == 8
        (V(a, 0, 0), V(0, b, 0), V(0, 0, c))
    elseif ibrav == 9
        (V(a / 2, b / 2, 0), V(-a / 2, b / 2, 0), V(0, 0, c))
    elseif ibrav == 10
        (V(a / 2, 0, c / 2), V(a / 2, b / 2, 0), V(0, b / 2, c / 2))
    elseif ibrav == 11
        (V(a / 2, b / 2, c / 2), V(-a / 2, b / 2, c / 2), V(-a / 2, -b / 2, c / 2))
    elseif ibrav == 12
        cγ = celldm[4]
        (V(a, 0, 0), V(b * cγ, b * sqrt(1 - cγ^2), 0), V(0, 0, c))
    elseif ibrav == -12
        cβ = celldm[5]
        (V(a, 0, 0), V(0, b, 0), V(c * cβ, 0, c * sqrt(1 - cβ^2)))
    elseif ibrav == 14
        cα, cβ, cγ = celldm[4], celldm[5], celldm[6]
        sγ = sqrt(1 - cγ^2)
        (V(a, 0, 0), V(b * cγ, b * sγ, 0),
         V(c * cβ, c * (cα - cβ * cγ) / sγ, c * sqrt(1 + 2cα * cβ * cγ - cα^2 - cβ^2 - cγ^2) / sγ))
    else
        throw(ArgumentError("ibrav = $ibrav is not supported (use 0–12, -3, -12 or 14)"))
    end
    return Lattice(vs...)
end

"""celldm(1..6) from the namelist, or from A, B, C, cosAB, cosAC, cosBC (Å)."""
function qe_celldm(vars, ibrav)
    cd = zeros(6)
    if haskey(vars, "a")
        A = qe_float(vars["a"])
        cd[1] = A / BOHR
        cd[2] = qe_float(get(vars, "b", vars["a"])) / A
        cd[3] = qe_float(get(vars, "c", vars["a"])) / A
        cab, cac, cbc = (qe_float(get(vars, k, "0")) for k in ("cosab", "cosac", "cosbc"))
        if ibrav == 14
            cd[4:6] .= (cbc, cac, cab)
        elseif ibrav == -12
            cd[5] = cac
        else
            cd[4] = cab
        end
    else
        for k in 1:6
            haskey(vars, "celldm($k)") && (cd[k] = qe_float(vars["celldm($k)"]))
        end
    end
    return cd
end

"""Card lines: the option of card `name` and its following lines (up to the next card)."""
function qe_card(lines, name)
    k = findfirst(l -> uppercase(split(l)[1]) == name, lines)
    k === nothing && return nothing, String[]
    opt = lowercase(strip(replace(lines[k][length(name)+1:end], r"[(){}=]" => " ")))
    body = String[]
    cards = ("ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS", "OCCUPATIONS",
             "CONSTRAINTS", "ATOMIC_FORCES", "ADDITIONAL_K_POINTS", "SOLVENTS", "HUBBARD", "&")
    for l in lines[k+1:end]
        any(c -> startswith(uppercase(l), c), cards) && break
        push!(body, l)
    end
    return opt, body
end

"""
    read_qe(path or io) -> Crystal

pw.x input: ibrav 0 with CELL_PARAMETERS (alat, bohr or angstrom) or ibrav
1–12, -3, -12, 14 with celldm or A, B, C, cos...; ATOMIC_POSITIONS in crystal,
alat, bohr or angstrom. Species labels map to elements via their name
(Fe1 → Fe).
"""
read_qe(path::AbstractString) = open(io -> read_qe(io; title=splitext(basename(path))[1]), path)

function read_qe(io::IO; title::AbstractString="")
    text = read(io, String)
    occursin(r"&system"i, text) || throw(ArgumentError("not a pw.x input (no &SYSTEM namelist)"))
    vars = qe_namelists(text)
    # Cards come after the last namelist.
    ends = collect(eachmatch(r"\n\s*/", text))
    cards = isempty(ends) ? text : text[ends[end].offset+length(ends[end].match):end]
    lines = filter(!isempty, [strip(replace(l, r"[!#].*" => "")) for l in split(cards, '\n')])
    ibrav = parse(Int, get(vars, "ibrav", "0"))
    cd = qe_celldm(vars, ibrav)
    alat = cd[1] * BOHR
    L = if ibrav == 0
        opt, body = qe_card(lines, "CELL_PARAMETERS")
        length(body) >= 3 || throw(ArgumentError("ibrav = 0 needs CELL_PARAMETERS"))
        # No unit: alat when celldm(1) / A is given, else Bohr (old pw.x convention).
        startswith(opt, "alat") && alat == 0 && throw(ArgumentError("CELL_PARAMETERS alat needs celldm(1) or A"))
        s = startswith(opt, "ang") ? 1.0 : startswith(opt, "bohr") ? BOHR : alat > 0 ? alat : BOHR
        Lattice((s * Vec3(qe_float.(split(body[i])[1:3])) for i in 1:3)...)
    else
        qe_lattice(ibrav, cd)
    end
    opt, body = qe_card(lines, "ATOMIC_POSITIONS")
    nat = parse(Int, get(vars, "nat", string(length(body))))
    length(body) >= nat || throw(ArgumentError("ATOMIC_POSITIONS lists $(length(body)) atoms, nat = $nat"))
    species = [species_token(split(l)[1]) for l in body[1:nat]]
    xyz = [Vec3(qe_float.(split(l)[2:4])) for l in body[1:nat]]
    isempty(opt) && (opt = "alat")
    startswith(opt, "crystal_sg") && throw(ArgumentError("ATOMIC_POSITIONS crystal_sg is not supported"))
    # if_pos flags (0 = fixed) after the coordinates.
    fixed = map(body[1:nat]) do l
        w = split(l)
        length(w) >= 7 ? Tuple(strip(t) == "0" for t in w[5:7]) : FREE_ATOM
    end
    startswith(opt, "crystal") && return Crystal(L, species, xyz; title, fixed)
    s = startswith(opt, "ang") ? 1.0 : startswith(opt, "bohr") ? BOHR : alat
    return with_constraints(crystal_from_cartesian(L, species, [s * r for r in xyz]; title); fixed)
end

"""
    read_qe_output(path or io) -> Crystal

pw.x output: the last CELL_PARAMETERS / ATOMIC_POSITIONS printed during a
relax or vc-relax (the final structure), otherwise the input structure
listed in the header.
"""
read_qe_output(path::AbstractString) = open(io -> read_qe_output(io; title=splitext(basename(path))[1]), path)

function read_qe_output(io::IO; title::AbstractString="")
    lines = readlines(io)
    any(l -> occursin("Program PWSCF", l), lines) || throw(ArgumentError("not a pw.x output (no 'Program PWSCF')"))
    num(s) = parse(Float64, s)
    alat = NaN
    L, species, pos, unit = nothing, Symbol[], Vec3[], "alat"
    k = 1
    while k <= length(lines)
        l = lines[k]
        if (m = match(r"lattice parameter \(alat\)\s*=\s*([\d.Ee+-]+)", l)) !== nothing
            alat = num(m[1]) * BOHR
        elseif occursin("crystal axes: (cart. coord. in units of alat)", l)
            L = Lattice((alat * Vec3(num.(split(match(r"\(([^)]*)\)\s*$", lines[k+i])[1]))) for i in 1:3)...)
            k += 3
        elseif occursin("positions (alat units)", l)
            empty!(species); empty!(pos)
            while k + 1 <= length(lines) && (m = match(r"^\s*\d+\s+(\S+)\s+tau\(\s*\d+\)\s*=\s*\(([^)]*)\)", lines[k+1])) !== nothing
                push!(species, species_token(m[1])); push!(pos, Vec3(num.(split(m[2]))))
                k += 1
            end
            unit = "alat"
        elseif startswith(strip(l), "CELL_PARAMETERS")
            s = (m = match(r"alat\s*=\s*([\d.Ee+-]+)", l)) !== nothing ? num(m[1]) * BOHR :
                occursin("angstrom", l) ? 1.0 : occursin("bohr", l) ? BOHR : alat
            L = Lattice((s * Vec3(num.(split(lines[k+i])[1:3])) for i in 1:3)...)
            k += 3
        elseif startswith(strip(l), "ATOMIC_POSITIONS")
            m = match(r"ATOMIC_POSITIONS\s*[({]?\s*(\w+)?", l)
            unit = m === nothing || m[1] === nothing ? "alat" : lowercase(m[1])
            empty!(species); empty!(pos)
            while k + 1 <= length(lines) && length(split(lines[k+1])) >= 4 && tryparse(Float64, split(lines[k+1])[2]) !== nothing
                w = split(lines[k+1])
                push!(species, species_token(w[1])); push!(pos, Vec3(num.(w[2:4])))
                k += 1
            end
        end
        k += 1
    end
    (L === nothing || isempty(pos)) && throw(ArgumentError("no cell or atoms found"))
    unit == "crystal" && return Crystal(L, species, pos; title)
    s = unit == "angstrom" ? 1.0 : unit == "bohr" ? BOHR : alat
    return crystal_from_cartesian(L, species, [s * r for r in pos]; title)
end
