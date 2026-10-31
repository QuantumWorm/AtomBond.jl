# Input files for ab-initio codes, plus a P1 CIF writer.
#
# k-point meshes come from a spacing in Å⁻¹ using the VASP KSPACING
# convention (reciprocal vectors include 2π): n_i = max(1, ceil(|b_i| / Δk)).
# Non-periodic directions and vacuum directions (slabs) get a single k-point.

"""
    vacuum_axes(c; min_gap=5.0) -> NTuple{3,Bool}

True along each lattice direction whose largest empty gap between atomic
planes (measured perpendicular to the planes, across the periodic boundary)
exceeds `min_gap` Å, e.g. the c axis of a slab.
"""
function vacuum_axes(c::Crystal; min_gap::Real=5.0)
    natoms(c) == 0 && return (true, true, true)
    d = plane_spacings(c.lattice)
    ntuple(3) do k
        z = sort(mod.([f[k] for f in c.frac], 1.0))
        gaps = [diff(z); z[1] + 1 - z[end]]
        maximum(gaps) * d[k] > min_gap
    end
end

"""
    kpoint_grid(c; kspacing=0.25, min_gap=5.0) -> (n1, n2, n3)

Monkhorst–Pack / Γ-centred mesh for a k-point spacing in Å⁻¹ (2π included,
as VASP's KSPACING; 0.25 is a reasonable default for semiconductors).
Directions that are not periodic or that contain more than `min_gap` Å of
vacuum get one k-point.
"""
function kpoint_grid(c::Crystal; kspacing::Real=0.25, min_gap::Real=5.0)
    kspacing > 0 || throw(ArgumentError("kspacing must be positive"))
    B = 2π * inv(c.lattice.matrix)'      # columns: reciprocal vectors
    vac = vacuum_axes(c; min_gap)
    ntuple(k -> c.pbc[k] && !vac[k] ? max(1, ceil(Int, norm(B[:, k]) / kspacing - 1e-8)) : 1, 3)
end

# Species in first-appearance order and the atoms grouped by species, the
# order VASP needs (and harmless for the other codes).
grouped_order(c::Crystal) = [i for s in species_list(c) for i in eachindex(c.species) if c.species[i] == s]

# ---------------------------------------------------------------------------
# VASP

"""
    write_kpoints(path, c; kspacing=0.25)

VASP KPOINTS file with a Γ-centred automatic mesh (Γ-centring is safe for
hexagonal cells as well).
"""
write_kpoints(path::AbstractString, c::Crystal; kw...) = open(io -> write_kpoints(io, c; kw...), path, "w")

function write_kpoints(io::IO, c::Crystal; kspacing::Real=0.25)
    n = kpoint_grid(c; kspacing)
    println(io, "Automatic mesh, k-spacing $(kspacing) 1/Å")
    println(io, "0")
    println(io, "Gamma")
    println(io, "  ", join(n, " "))
    println(io, "  0 0 0")
end

"""
    write_vasp_inputs(dir, c; kspacing=0.25) -> Vector of paths

POSCAR and KPOINTS in `dir` (created if needed). POTCAR is licensed and must
be assembled by the user in the species order of the POSCAR, which is also
written to POTCAR.order for reference.
"""
function write_vasp_inputs(dir::AbstractString, c::Crystal; kspacing::Real=0.25)
    mkpath(dir)
    paths = [joinpath(dir, "POSCAR"), joinpath(dir, "KPOINTS"), joinpath(dir, "POTCAR.order")]
    write_poscar(paths[1], c)
    write_kpoints(paths[2], c; kspacing)
    write(paths[3], join(string.(species_list(c)), " ") * "\n")
    if has_magmoms(c)
        # Initial moments in POSCAR order (unknown = 0), as runs "n*value".
        push!(paths, joinpath(dir, "INCAR.spin"))
        m = [isnan(c.magmoms[i]) ? 0.0 : c.magmoms[i] for s in species_list(c) for i in eachindex(c.species)
             if c.species[i] == s]
        runs = Tuple{Int,Float64}[]
        for x in m
            !isempty(runs) && runs[end][2] == x ? (runs[end] = (runs[end][1] + 1, x)) : push!(runs, (1, x))
        end
        open(paths[end], "w") do io
            println(io, "# Spin polarization with the initial magnetic moments (μB) of AtomBond; merge into INCAR")
            println(io, "ISPIN = 2")
            println(io, "MAGMOM = ", join((n == 1 ? @sprintf("%g", x) : @sprintf("%d*%g", n, x) for (n, x) in runs), " "))
        end
    end
    if has_constraints(c)
        # Fixed atoms are in the POSCAR (selective dynamics); the cell needs ISIF.
        push!(paths, joinpath(dir, "INCAR.relax"))
        isif, note = cell_constraints_vasp(c)
        open(paths[end], "w") do io
            println(io, "# Geometry optimization settings for the constraints (AtomBond); merge into INCAR")
            println(io, "# Cell: ", describe_cell_constraints(c.cell_constraints), ". ", note)
            println(io, "IBRION = 2\nNSW = 100\nEDIFFG = -0.01")
            println(io, isif === nothing ? "ISIF = 3   # not exact, see above" : "ISIF = $isif")
        end
    end
    return paths
end

# ---------------------------------------------------------------------------
# Quantum ESPRESSO

"""
    write_qe(path, c; kspacing=0.25, ecutwfc=40, pseudo="{El}.UPF", calculation=automatic)

pw.x input with `ibrav = 0`, the cell in Å and positions in crystal
coordinates. `pseudo` is a file-name pattern in which `{El}` becomes the
element symbol. `ecutwfc` is in Ry; check it against your pseudopotentials.

Constraints: fixed atoms get `if_pos` flags (0 = fixed; pw.x applies them to
Cartesian force components, which equals fixing fractional coordinates only
when the cell axes are along x, y, z; a fully fixed atom is always exact).
`calculation` is by default 'scf' without constraints, 'relax' for a fixed
cell and 'vc-relax' with `cell_dofree` otherwise (see `cell_constraints_qe`).
"""
write_qe(path::AbstractString, c::Crystal; kw...) = open(io -> write_qe(io, c; kw...), path, "w")

function write_qe(io::IO, c::Crystal; kspacing::Real=0.25, ecutwfc::Real=40,
                  pseudo::AbstractString="{El}.UPF", calculation::Union{Nothing,AbstractString}=nothing)
    spc = species_list(c)
    n = kpoint_grid(c; kspacing)
    prefix = formula(c)
    auto, dofree, note = has_constraints(c) ? cell_constraints_qe(c) : ("scf", nothing, "")
    calculation = something(calculation, auto)
    println(io, "! ", isempty(c.title) ? prefix : c.title, " — written by AtomBond")
    isempty(note) || println(io, "! Constraints: ", describe_cell_constraints(c.cell_constraints), "; ", note)
    # if_pos acts on Cartesian force components: fixing fractional coordinate k
    # is exact when the reciprocal direction k lies along Cartesian axis k.
    B = inv(c.lattice.matrix)
    exact(k) = (n = normalize(Vec3(B[k, :])); abs(abs(n[k]) - 1) < 1e-8)
    any(f -> any(f) && !all(f) && any(k -> f[k] && !exact(k), 1:3), c.fixed) &&
        println(io, "! Note: if_pos fixes Cartesian force components; some partly fixed atoms are approximate in this cell")
    println(io, "&CONTROL\n  calculation = '$calculation'\n  prefix = '$prefix'\n",
                "  pseudo_dir = './'\n  outdir = './tmp'\n/")
    println(io, "&SYSTEM\n  ibrav = 0\n  nat = $(natoms(c))\n  ntyp = $(length(spc))\n",
                "  ecutwfc = $(Float64(ecutwfc))\n/")
    println(io, "&ELECTRONS\n  conv_thr = 1.0d-8\n/")
    calculation in ("relax", "vc-relax", "md", "vc-md") && println(io, "&IONS\n/")
    if calculation in ("vc-relax", "vc-md")
        println(io, "&CELL\n", dofree === nothing ? "" : "  cell_dofree = '$dofree'\n", "/")
    end
    println(io, "ATOMIC_SPECIES")
    for s in spc
        @printf(io, "  %-3s %10.4f  %s\n", s, atomic_mass(s), replace(pseudo, "{El}" => string(s)))
    end
    println(io, "CELL_PARAMETERS angstrom")
    for v in cell_vectors(c.lattice)
        @printf(io, "  %16.10f %16.10f %16.10f\n", v...)
    end
    println(io, "ATOMIC_POSITIONS crystal")
    flags = any(any, c.fixed)
    for i in grouped_order(c)
        @printf(io, "  %-3s %16.10f %16.10f %16.10f", c.species[i], c.frac[i]...)
        flags && print(io, "  ", join((x ? 0 : 1 for x in c.fixed[i]), " "))
        println(io)
    end
    println(io, "K_POINTS automatic")
    println(io, "  ", join(n, " "), " 0 0 0")
end

# ---------------------------------------------------------------------------
# CASTEP

"""
    write_castep(path, c; kspacing=0.25)

CASTEP `.cell` file: Cartesian lattice in Å, fractional positions and an
MP grid (CASTEP generates pseudopotentials on the fly when no SPECIES_POT
block is given). Write the matching `.param` file yourself.
"""
write_castep(path::AbstractString, c::Crystal; kw...) = open(io -> write_castep(io, c; kw...), path, "w")

function write_castep(io::IO, c::Crystal; kspacing::Real=0.25)
    println(io, "! ", isempty(c.title) ? formula(c) : c.title, " — written by AtomBond\n")
    println(io, "%BLOCK LATTICE_CART\nang")
    for v in cell_vectors(c.lattice)
        @printf(io, "  %16.10f %16.10f %16.10f\n", v...)
    end
    println(io, "%ENDBLOCK LATTICE_CART\n")
    println(io, "%BLOCK POSITIONS_FRAC")
    for i in grouped_order(c)
        @printf(io, "  %-3s %16.10f %16.10f %16.10f", c.species[i], c.frac[i]...)
        isnan(c.magmoms[i]) || @printf(io, "  SPIN=%.4f", c.magmoms[i])      # initial moment, μB
        println(io)
    end
    println(io, "%ENDBLOCK POSITIONS_FRAC\n")
    write_castep_constraints(io, c)
    println(io, "KPOINTS_MP_GRID ", join(kpoint_grid(c; kspacing), " "))
    println(io, "SYMMETRY_GENERATE")
end

"""
CELL_CONSTRAINTS (a b c / α β γ) and IONIC_CONSTRAINTS. A fixed fractional
coordinate k of an atom is the linear constraint r · n = const with n along
the reciprocal vector k (Cartesian, as CASTEP wants); ions are numbered
within their species in the order of POSITIONS_FRAC.
"""
function write_castep_constraints(io::IO, c::Crystal)
    if c.cell_constraints != FREE_CELL
        cc = c.cell_constraints
        println(io, "! ", describe_cell_constraints(cc))
        println(io, "%BLOCK CELL_CONSTRAINTS\n  ", join(cc[1:3], " "), "\n  ", join(cc[4:6], " "),
                "\n%ENDBLOCK CELL_CONSTRAINTS\n")
    end
    any(any, c.fixed) || return
    B = inv(c.lattice.matrix)                      # rows: reciprocal directions of a, b, c
    println(io, "%BLOCK IONIC_CONSTRAINTS")
    n, count = 0, Dict{Symbol,Int}()
    for i in grouped_order(c)
        s = c.species[i]
        k = count[s] = get(count, s, 0) + 1
        for axis in 1:3
            c.fixed[i][axis] || continue
            v = normalize(Vec3(B[axis, :]))
            n += 1
            @printf(io, "  %3d  %-3s %4d  %14.10f %14.10f %14.10f\n", n, s, k, map(x -> abs(x) < 1e-12 ? 0.0 : x, v)...)
        end
    end
    println(io, "%ENDBLOCK IONIC_CONSTRAINTS\n")
end

# ---------------------------------------------------------------------------
# CIF

"""
    write_cif(path, c)

P1 CIF (all atoms listed, no symmetry operations beyond the identity).
"""
write_cif(path::AbstractString, c::Crystal) = open(io -> write_cif(io, c), path, "w")

function write_cif(io::IO, c::Crystal)
    name = replace(formula(c), r"\W" => "_")
    a, b, cc, α, β, γ = lattice_parameters(c.lattice)
    println(io, "# written by AtomBond")
    println(io, "data_", name)
    isempty(c.title) || println(io, "_chemical_name_systematic  '", replace(c.title, "'" => ""), "'")
    @printf(io, "_cell_length_a     %.6f\n_cell_length_b     %.6f\n_cell_length_c     %.6f\n", a, b, cc)
    @printf(io, "_cell_angle_alpha  %.6f\n_cell_angle_beta   %.6f\n_cell_angle_gamma  %.6f\n", α, β, γ)
    println(io, "_symmetry_space_group_name_H-M  'P 1'\n_symmetry_Int_Tables_number  1")
    println(io, "loop_\n_symmetry_equiv_pos_as_xyz\n  'x, y, z'")
    println(io, "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z")
    for i in eachindex(c.species)
        @printf(io, "  %-6s %-3s %12.8f %12.8f %12.8f\n", atom_label(c, i), c.species[i], c.frac[i]...)
    end
end

# ---------------------------------------------------------------------------
# One entry point for the GUI

const EXPORT_FORMATS = [
    ("VASP (POSCAR + KPOINTS)", :vasp),
    ("Quantum ESPRESSO (.in)", :qe),
    ("CASTEP (.cell)", :castep),
    ("CIF (P1)", :cif),
    ("Extended XYZ", :xyz),
]

"""
    export_inputs(c, format, dir, name; kspacing=0.25, ecutwfc=40, pseudo="{El}.UPF") -> paths

Write `c` for `format` ∈ (:vasp, :qe, :castep, :cif, :xyz) into `dir`. VASP
files go into a sub-folder `name`; the others are `dir/name.<ext>`.
"""
function export_inputs(c::Crystal, format::Symbol, dir::AbstractString, name::AbstractString;
                       kspacing::Real=0.25, ecutwfc::Real=40, pseudo::AbstractString="{El}.UPF")
    isempty(strip(name)) && throw(ArgumentError("empty file name"))
    mkpath(dir)
    format === :vasp && return write_vasp_inputs(joinpath(dir, name), c; kspacing)
    ext = Dict(:qe => ".in", :castep => ".cell", :cif => ".cif", :xyz => ".xyz")
    haskey(ext, format) || throw(ArgumentError("unknown export format $format"))
    path = joinpath(dir, name * ext[format])
    format === :qe ? write_qe(path, c; kspacing, ecutwfc, pseudo) :
    format === :castep ? write_castep(path, c; kspacing) :
    format === :cif ? write_cif(path, c) : write_xyz(path, c)
    return [path]
end
