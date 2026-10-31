# Pick a reader or writer from the file name.

"""
    structure_format(path) -> Symbol or nothing

Which reader [`load_structure`](@ref) would use for this file name, or
`nothing`: `:cif`, `:xyz`, `:poscar`, `:xdatcar`, `:vasprun`, `:castep_cell`,
`:castep_geom`, `:castep_out`, `:qe_in`, `:qe_out`, `:aims`, `:siesta`, `:abinit`, `:xsf`,
`:cube`, `:pdb`, `:lammps`, `:car`, `:xsd`, `:res`.
"""
function structure_format(path::AbstractString)
    base = lowercase(basename(path))
    ext = lowercase(splitext(base)[2])
    startswith(base, "geometry.in") && return :aims
    occursin("xdatcar", base) && return :xdatcar
    startswith(base, "vasprun") && ext == ".xml" && return :vasprun
    (ext in (".vasp", ".poscar") || startswith(base, "poscar") || startswith(base, "contcar")) && return :poscar
    startswith(base, "data.") && return :lammps
    return get(STRUCTURE_EXTENSIONS, ext, nothing)
end

const STRUCTURE_EXTENSIONS = Dict(
    ".cif" => :cif, ".xyz" => :xyz, ".extxyz" => :xyz, ".cell" => :castep_cell, ".geom" => :castep_geom,
    ".castep" => :castep_out,
    ".in" => :qe_in, ".pwi" => :qe_in, ".pwo" => :qe_out, ".out" => :qe_out,
    ".fdf" => :siesta, ".abi" => :abinit, ".abinit" => :abinit,
    ".xsf" => :xsf, ".cube" => :cube, ".cub" => :cube, ".pdb" => :pdb, ".ent" => :pdb,
    ".lmp" => :lammps, ".lammps" => :lammps, ".data" => :lammps,
    ".car" => :car, ".xsd" => :xsd, ".res" => :res, ".ins" => :res)

"""Reader for a `.in` file: pw.x, FHI-aims or ABINIT, by content."""
function read_in_file(path)
    text = read(path, String)
    occursin(r"&system"i, text) && return read_qe(path)
    occursin(r"^\s*(lattice_vector|atom_frac|atom)\s"m, text) && return read_aims(path)
    occursin(r"^\s*(acell|xred|xcart|xangst|typat|znucl)\b"m, text) && return read_abinit(path)
    throw(ArgumentError("$(basename(path)): not a pw.x, FHI-aims or ABINIT input"))
end

"""
    load_structure(path) -> Crystal

Reads (format from the file name, see [`structure_format`](@ref)):

- VASP: POSCAR/CONTCAR/`.vasp`, XDATCAR (last frame), vasprun.xml (final structure)
- CASTEP: `.cell`, `.geom` (last frame, i.e. the relaxed structure),
  `.castep` output (final structure with Mulliken charges)
- Quantum ESPRESSO: pw.x input `.in`/`.pwi`, pw.x output `.out`/`.pwo` (last structure)
- FHI-aims `geometry.in`, SIESTA `.fdf`, ABINIT `.abi`/`.abinit` (and `.in` by content)
- CIF, XYZ / extended XYZ, XSF (XCrySDen, VESTA), Gaussian cube, PDB
- LAMMPS data (`.lmp`, `.data`, `data.*`), Materials Studio `.car` / `.xsd`, SHELX/AIRSS `.res`/`.ins`
"""
function load_structure(path::AbstractString)
    isfile(path) || throw(ArgumentError("no such file: $path"))
    fmt = structure_format(path)
    fmt === :cif && return read_cif(path)
    fmt === :xyz && return read_xyz(path)
    fmt === :poscar && return read_poscar(path)
    fmt === :xdatcar && return read_xdatcar(path)
    fmt === :vasprun && return read_vasprun(path)
    fmt === :castep_cell && return read_castep_cell(path)
    fmt === :castep_geom && return read_castep_geom(path)
    fmt === :castep_out && return read_castep_output(path)
    fmt === :qe_in && return read_in_file(path)
    fmt === :qe_out && return read_qe_output(path)
    fmt === :aims && return read_aims(path)
    fmt === :siesta && return read_siesta(path)
    fmt === :abinit && return read_abinit(path)
    fmt === :xsf && return read_xsf(path)
    fmt === :cube && return read_cube(path)
    fmt === :pdb && return read_pdb(path)
    fmt === :lammps && return read_lammps_data(path)
    fmt === :car && return read_car(path)
    fmt === :xsd && return read_xsd(path)
    fmt === :res && return read_res(path)
    throw(ArgumentError("don't know how to read '$(basename(path))' (see `?load_structure` for the formats)"))
end

"""
    save_structure(path, crystal; kwargs...)

Writer chosen by the file name: `.cif` (P1), POSCAR/CONTCAR/`.vasp`,
`.xyz`/`.extxyz`, `.in`/`.pwi` (Quantum ESPRESSO, see [`write_qe`](@ref)),
`.cell` (CASTEP, see [`write_castep`](@ref)), `geometry.in` (FHI-aims),
`.xsf`, `.pdb`, `.lmp`/`.data`/`.lammps` (LAMMPS data), `.car` (Materials
Studio), `.res` (SHELX/AIRSS). Keyword arguments (e.g. `kspacing`) go to the
QE / CASTEP writers.
"""
function save_structure(path::AbstractString, c::Crystal; kwargs...)
    base = lowercase(basename(path))
    ext = lowercase(splitext(path)[2])
    if startswith(base, "geometry.in")
        write_aims(path, c)
    elseif ext == ".cif"
        write_cif(path, c)
    elseif ext in (".in", ".pwi")
        write_qe(path, c; kwargs...)
    elseif ext == ".cell"
        write_castep(path, c; kwargs...)
    elseif ext in (".xyz", ".extxyz")
        write_xyz(path, c)
    elseif ext in (".vasp", ".poscar") || startswith(base, "poscar") || startswith(base, "contcar")
        write_poscar(path, c)
    elseif ext == ".xsf"
        write_xsf(path, c)
    elseif ext in (".pdb", ".ent")
        write_pdb(path, c)
    elseif ext in (".lmp", ".data", ".lammps") || startswith(base, "data.")
        write_lammps_data(path, c)
    elseif ext == ".car"
        write_car(path, c)
    elseif ext in (".res", ".ins")
        write_res(path, c)
    else
        throw(ArgumentError("don't know how to write '$(basename(path))'"))
    end
    return path
end
