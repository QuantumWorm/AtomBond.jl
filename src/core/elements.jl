# Element data, read from data/elements.csv (one row per element, Z 1–103).
#
# The table holds the published values with their sources (see its header and
# scripts/build_element_table.jl):
#   mass                 NIST / IUPAC standard atomic weights
#   covalent_radius      Cordero et al. 2008      → bonds, ball size
#   vdw_radius_bondi     Bondi 1964 (+ Mantina)   → space-filling (CPK) size
#   vdw_radius_alvarez   Alvarez 2013             → CPK size where Bondi has none
#   color_jmol           Jmol colors (O red, N blue, C grey, H white, ...)
#
# `Element` holds the values the program uses, picked from those columns.
# Any column (also new ones, e.g. electronegativity) can be read with
# `element_property`. Values can be overridden without touching the program:
# ~/.atombond/elements.csv (or $ATOMBOND_HOME/elements.csv) has the same
# format; its rows replace or add columns of the named elements, and the extra
# columns `color` (hex) and `vdw_radius` take precedence over the defaults.
# `use_vesta_colors!` takes colors from a VESTA elements.ini.

struct Element
    symbol::Symbol
    Z::Int
    name::String
    mass::Float64              # u (standard atomic weight)
    covalent_radius::Float64   # Å
    vdw_radius::Float64        # Å
    color::RGBf
end

Base.show(io::IO, e::Element) = print(io, "Element(", e.symbol, ", Z=", e.Z, ")")

hexcolor(h::AbstractString) = RGBf(parse(UInt8, h[1:2]; base=16) / 255,
                                   parse(UInt8, h[3:4]; base=16) / 255,
                                   parse(UInt8, h[5:6]; base=16) / 255)

const ELEMENT_FILE = normpath(joinpath(@__DIR__, "..", "..", "data", "elements.csv"))
include_dependency(ELEMENT_FILE)

"""Rows of an element table (CSV with a header line; # starts a comment line) as column => text."""
function read_element_table(path::AbstractString)
    lines = filter(l -> !isempty(strip(l)) && !startswith(lstrip(l), '#'), readlines(path))
    isempty(lines) && return Dict{Symbol,String}[]
    header = Symbol.(strip.(split(lines[1], ',')))
    :symbol in header || throw(ArgumentError("$path: the header needs a 'symbol' column"))
    return map(lines[2:end]) do l
        w = strip.(split(l, ','))
        length(w) == length(header) || throw(ArgumentError("$path: '$l' has $(length(w)) values for $(length(header)) columns"))
        Dict(zip(header, String.(w)))
    end
end

# symbol => column => text (the built-in table, plus overrides).
const ELEMENT_TABLE = Dict{Symbol,Dict{Symbol,String}}(Symbol(r[:symbol]) => r for r in read_element_table(ELEMENT_FILE))

"""
    element_property(sym, column) -> Float64, String or nothing

Any column of the element table for an element: `:mass`, `:covalent_radius`,
`:vdw_radius_bondi`, `:vdw_radius_alvarez`, `:color_jmol`, `:name`, or a
column added in an override file. Numbers are returned as Float64, empty
values as `nothing`.
"""
function element_property(s, column::Symbol)
    row = get(ELEMENT_TABLE, element(s).symbol, nothing)
    x = row === nothing ? "" : get(row, column, "")
    isempty(x) && return nothing
    v = tryparse(Float64, x)
    return v === nothing ? x : v
end

const DEFAULT_COVALENT_RADIUS = 1.6     # Å, for elements without a tabulated value (Z > 96)

function make_element(row::Dict{Symbol,String})
    num(k) = (x = get(row, k, ""); isempty(x) ? nothing : parse(Float64, x))
    text(k) = (x = get(row, k, ""); isempty(x) ? nothing : x)
    r = something(num(:covalent_radius), DEFAULT_COVALENT_RADIUS)
    Element(Symbol(row[:symbol]), parse(Int, row[:Z]), something(text(:name), row[:symbol]),
            something(num(:mass), NaN), r,
            something(num(:vdw_radius), num(:vdw_radius_bondi), num(:vdw_radius_alvarez), r + 0.8),
            hexcolor(something(text(:color), text(:color_jmol), "FF1493")))
end

const ELEMENTS = Dict{Symbol,Element}()
const ELEMENTS_BY_Z = Dict{Int,Element}()

"""Rebuild the element lookup from `ELEMENT_TABLE` (after overrides)."""
function refresh_elements!()
    empty!(ELEMENTS); empty!(ELEMENTS_BY_Z)
    for row in values(ELEMENT_TABLE)
        e = make_element(row)
        ELEMENTS[e.symbol] = e
        ELEMENTS_BY_Z[e.Z] = e
    end
end
refresh_elements!()

"""
    override_elements!(rows) / override_elements!(path)

Merge rows (column => text, with a `symbol` column) into the element table:
given columns replace the built-in ones; a new symbol needs `Z`. Affects
structures drawn afterwards.
"""
function override_elements!(rows::AbstractVector)
    for r in rows
        sym = Symbol(r[:symbol])
        base = get(ELEMENT_TABLE, sym, Dict{Symbol,String}())
        ELEMENT_TABLE[sym] = merge(base, Dict(k => v for (k, v) in r if !isempty(v)))
    end
    refresh_elements!()
end
override_elements!(path::AbstractString) = override_elements!(read_element_table(path))

"""
    use_vesta_colors!(path; radii=false)

Take element colors (and, with `radii=true`, the van der Waals radii used
for CPK) from a VESTA `elements.ini` (lines: number, symbol, atomic radius,
van der Waals radius, ionic radius, R G B in 0–1).
"""
function use_vesta_colors!(path::AbstractString; radii::Bool=false)
    rows = Dict{Symbol,String}[]
    for l in eachline(path)
        w = split(l)
        (length(w) >= 8 && haskey(ELEMENT_TABLE, Symbol(w[2]))) || continue      # skips D, XX, ...
        rgb = round.(Int, 255 .* clamp.(parse.(Float64, w[6:8]), 0, 1))
        r = Dict(:symbol => String(w[2]), :color => join(string.(rgb; base=16, pad=2)))
        radii && (r[:vdw_radius] = String(w[4]))
        push!(rows, r)
    end
    isempty(rows) && throw(ArgumentError("$path: no element lines (not a VESTA elements.ini?)"))
    override_elements!(rows)
    return length(rows)
end

"""Apply ~/.atombond/elements.csv (or \$ATOMBOND_HOME/elements.csv) if it exists."""
function load_user_elements!()
    path = joinpath(get(ENV, "ATOMBOND_HOME", joinpath(homedir(), ".atombond")), "elements.csv")
    isfile(path) || return nothing
    try
        override_elements!(path)
    catch err
        @warn "Ignoring $path" exception = err
    end
    return path
end

# Fallback for unknown species (e.g. dummy atoms "X"): magenta, medium size.
const UNKNOWN_ELEMENT = Element(:X, 0, "unknown", NaN, 1.0, 1.8, RGBf(1.0, 0.08, 0.58))

"""
    element(sym) -> Element

Look up an element by symbol (`:Si`, `"si"`, `"Si1"`) or atomic number.
Unknown symbols return a magenta placeholder rather than throwing, so that
odd labels in input files still render.
"""
element(z::Integer) = get(ELEMENTS_BY_Z, z, UNKNOWN_ELEMENT)
element(s::Symbol) = get(ELEMENTS, s, UNKNOWN_ELEMENT)
element(s::AbstractString) = element(normalize_symbol(s))
element(e::Element) = e

"""Van der Waals radius (Å) used for CPK: Bondi, else Alvarez, else covalent + 0.8 Å."""
vdw_radius(s) = element(s).vdw_radius

"""Standard atomic weight (u); throws for unknown species."""
function atomic_mass(s)
    e = element(s)
    (e.Z == 0 || isnan(e.mass)) && throw(ArgumentError("no atomic mass for '$s'"))
    return e.mass
end

# Non-metals and metalloids; everything else counts as a metal for bonding.
const NONMETALS = Set([:H, :He, :B, :C, :N, :O, :F, :Ne, :Si, :P, :S, :Cl, :Ar,
                       :Ge, :As, :Se, :Br, :Kr, :Sb, :Te, :I, :Xe])

ismetal(s::Symbol) = s ∉ NONMETALS && s !== :X

"""
    normalize_symbol(s) -> Symbol

Turn labels like `"SI"`, `"si"`, `"Si1"`, `"O2-"` into a clean element symbol.
"""
function normalize_symbol(s::AbstractString)
    letters = filter(isletter, strip(s))
    isempty(letters) && return :X
    # Try a two-letter symbol first, then one letter.
    if length(letters) >= 2
        two = Symbol(uppercase(letters[1]) * lowercase(letters[2]))
        haskey(ELEMENTS, two) && return two
    end
    return Symbol(uppercase(letters[1:1]))
end
