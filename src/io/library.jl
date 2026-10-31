# Structure library: known crystals and molecules with their sources.
#
# The built-in library is data/library/: an index, library.toml, with one
# [[structure]] entry per structure (name, kind, category, formula, file,
# description, source, url, reference, doi, note), and the files it names
# (crystals: CIFs from the Crystallography Open Database; molecules: XYZ
# files of NIST CCCBDB experimental geometries). scripts/build_library.jl
# rebuilds it from those sources.
#
# Your own structures: a library.toml in ~/.atombond/library/ (or
# $ATOMBOND_HOME/library/) in the same format is added to the list, with
# files relative to that folder.

const LIBRARY_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "library"))

"""One structure of the library; `file` is an absolute path."""
struct LibraryEntry
    name::String
    kind::Symbol            # :crystal or :molecule
    category::String
    formula::String
    file::String
    description::String
    source::String
    url::String
    reference::String
    doi::String
    note::String
end

Base.show(io::IO, e::LibraryEntry) = print(io, "LibraryEntry(\"", e.name, "\", ", e.kind, ")")

"""Entries of one library.toml; `files` are resolved relative to its folder."""
function read_library_index(path::AbstractString)
    dir = dirname(path)
    entries = get(TOML.parsefile(path), "structure", Any[])
    return map(entries) do d
        field(k) = string(get(d, k, ""))
        haskey(d, "name") && haskey(d, "file") || throw(ArgumentError("$path: every [[structure]] needs a name and a file"))
        kind = Symbol(get(d, "kind", "crystal"))
        kind in (:crystal, :molecule) || throw(ArgumentError("$path: kind of '$(d["name"])' must be \"crystal\" or \"molecule\""))
        LibraryEntry(field("name"), kind, field("category"), field("formula"), normpath(joinpath(dir, d["file"])),
                     field("description"), field("source"), field("url"), field("reference"), field("doi"), field("note"))
    end
end

"""Folder of the user's own library (may not exist)."""
user_library_dir() = joinpath(get(ENV, "ATOMBOND_HOME", joinpath(homedir(), ".atombond")), "library")

"""
    structure_library(; kind=nothing) -> Vector{LibraryEntry}

All structures of the built-in library, followed by those of the user's
library (~/.atombond/library/library.toml), if any. `kind=:crystal` or
`:molecule` keeps one kind.
"""
function structure_library(; kind::Union{Nothing,Symbol}=nothing)
    entries = read_library_index(joinpath(LIBRARY_DIR, "library.toml"))
    user = joinpath(user_library_dir(), "library.toml")
    if isfile(user)
        try
            append!(entries, read_library_index(user))
        catch err
            @warn "Ignoring $user" exception = err
        end
    end
    return kind === nothing ? entries : filter(e -> e.kind === kind, entries)
end

"""The library entry called `name` (or with that formula), ignoring case."""
function library_entry(name::AbstractString)
    entries = structure_library()
    key = lowercase(strip(name))
    i = something(findfirst(e -> lowercase(e.name) == key, entries),
                   findfirst(e -> lowercase(e.formula) == key, entries), Some(nothing))
    i === nothing && throw(ArgumentError("no structure '$name' in the library (see structure_library())"))
    return entries[i]
end

"""
    load_library(name or entry) -> Crystal

A structure from the library, by name (`"Silicon"`) or formula (`"H2O"`).
Molecules come in a box with 5 Å of space around them, as XYZ files do.
"""
load_library(name::AbstractString) = load_library(library_entry(name))

function load_library(e::LibraryEntry)
    c = load_structure(e.file)
    return Crystal(c.lattice, c.species, c.frac; pbc=c.pbc, title=e.name, fixed=c.fixed, cell_constraints=c.cell_constraints,
                   charges=c.charges, magmoms=c.magmoms)
end
