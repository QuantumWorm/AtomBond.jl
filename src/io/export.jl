# Input files for ab-initio codes (Quantum ESPRESSO, ABINIT, SIESTA, VASP,
# CASTEP), plus a P1 CIF writer.
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
n = ⌈|b|/Δk⌉; 0.25 is a reasonable default for semiconductors).
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
# ABINIT

"""Whether fixing fractional coordinate `k` equals fixing Cartesian component `k` in this cell."""
fixed_axis_is_cartesian(c::Crystal, k) = (n = normalize(Vec3(inv(c.lattice.matrix)[k, :])); abs(abs(n[k]) - 1) < 1e-8)

"""
    abinit_optcell(c) -> (optcell, note)

ABINIT `optcell` for the cell constraints: 0 cell fixed, 1 volume only,
2 everything, 6 only the length of c (acell(3)); others are written as 2.
"""
function abinit_optcell(c::Crystal)
    cc = c.cell_constraints
    cc == (0, 0, 0, 0, 0, 0) && return 0, "optcell 0: cell fixed"
    cc == FREE_CELL && return 2, "optcell 2: cell and atoms relax"
    cc == (1, 1, 1, 0, 0, 0) && return 1, "optcell 1: volume only"
    cc == (0, 0, 1, 0, 0, 0) && return 6, "optcell 6: only the length of c"
    return 2, "ABINIT has no optcell for these constraints; written as 2 (everything relaxes, edit by hand)"
end

"""
    write_abinit(path, c; kspacing=0.25, ecutwfc=40, pseudo="{El}.psp8")

ABINIT input (`.abi`, ABINIT 9 or newer): the cell (`acell` 1 Å with `rprim`
rows as the lattice vectors), `natom`, `ntypat`, `znucl`, `typat`, reduced
coordinates `xred`, a Γ-centred `ngkpt` mesh, `ecut` in Ha (= `ecutwfc`/2,
as `ecutwfc` is in Ry like Quantum ESPRESSO's) and `pseudos` from the
pattern (`{El}` becomes the element symbol). Initial magnetic moments give
`nsppol 2` and `spinat`.

Constraints: a relaxation (`ionmov 2`) with `optcell` for the cell (see
`abinit_optcell`); fully fixed atoms are listed in `iatfix`, and a fixed
fractional coordinate along a Cartesian axis in `iatfixx`/`iatfixy`/`iatfixz`
(other partly fixed atoms are noted in a comment).
"""
write_abinit(path::AbstractString, c::Crystal; kw...) = open(io -> write_abinit(io, c; kw...), path, "w")

function write_abinit(io::IO, c::Crystal; kspacing::Real=0.25, ecutwfc::Real=40, pseudo::AbstractString="{El}.psp8")
    spc = species_list(c)
    title = isempty(c.title) ? formula(c) : c.title
    println(io, "# ", title, " — written by AtomBond")
    println(io, "\n# Cell (lattice vectors as the rows of rprim, in Å)\nacell 3*1.0 angstrom\nrprim")
    for v in cell_vectors(c.lattice)
        @printf(io, "  %16.10f %16.10f %16.10f\n", v...)
    end
    println(io, "\n# Atoms\nnatom ", natoms(c), "\nntypat ", length(spc))
    println(io, "znucl ", join((element(s).Z for s in spc), " "))
    println(io, "typat ", join((findfirst(==(s), spc) for s in c.species), " "))
    println(io, "xred")
    for f in c.frac
        @printf(io, "  %16.10f %16.10f %16.10f\n", f...)
    end
    if any(!isnan, c.magmoms)
        println(io, "\n# Initial magnetic moments (μB, along z)\nnsppol 2\nspinat")
        foreach(m -> @printf(io, "  0.0 0.0 %.4f\n", isnan(m) ? 0.0 : m), c.magmoms)
    end
    n = kpoint_grid(c; kspacing)
    println(io, "\n# k-points: Γ-centred mesh for a spacing of $(kspacing) 1/Å")
    println(io, "kptopt 1\nngkpt ", join(n, " "), "\nnshiftk 1\nshiftk 0.0 0.0 0.0")
    println(io, "\n# Basis and pseudopotentials (check ecut against your pseudopotentials)")
    @printf(io, "ecut %.2f   # Ha (= %.2f Ry)\n", ecutwfc / 2, ecutwfc)
    println(io, "pp_dirpath \"./\"\npseudos \"", join((replace(pseudo, "{El}" => string(s)) for s in spc), ", "), "\"")
    println(io, "\n# SCF\nnstep 100\ntoldfe 1.0d-10")
    has_constraints(c) || return
    optcell, note = abinit_optcell(c)
    println(io, "\n# Relaxation: ", describe_cell_constraints(c.cell_constraints), "; ", note)
    println(io, "ionmov 2\nntime 100\ntolmxf 5.0d-5")
    optcell == 0 || println(io, "optcell ", optcell, "\ndilatmx 1.05\necutsm 0.5")
    full = [i for i in eachindex(c.fixed) if all(c.fixed[i])]
    isempty(full) || println(io, "natfix ", length(full), "\niatfix ", join(full, " "))
    approx = false
    for (k, x) in enumerate(("x", "y", "z"))
        ids = [i for i in eachindex(c.fixed) if c.fixed[i][k] && !all(c.fixed[i])]
        isempty(ids) && continue
        if fixed_axis_is_cartesian(c, k)
            println(io, "natfix$x ", length(ids), "\niatfix$x ", join(ids, " "))
        else
            approx = true
        end
    end
    approx && println(io, "# Note: some partly fixed atoms are not written (their fixed fractional directions " *
                          "are not along x, y or z in this cell)")
end

# ---------------------------------------------------------------------------
# SIESTA

"""
    write_siesta(path, c; kspacing=0.25)

SIESTA `.fdf` input: the species (pseudopotentials are read from
`<Label>.psf` or `.psml`), the lattice vectors in Å, fractional coordinates,
a Monkhorst–Pack k-grid, a DZP basis, PBE and `MeshCutoff 300 Ry` (adjust
it). Initial magnetic moments give `Spin polarized` and `DM.InitSpin`.

Constraints: a conjugate-gradient relaxation (`MD.TypeOfRun CG`) with
`MD.VariableCell` when the cell may change (SIESTA relaxes the whole cell;
other cell constraints are noted in a comment), and fully fixed atoms in
`Geometry.Constraints` (partly fixed atoms are noted in a comment).
"""
write_siesta(path::AbstractString, c::Crystal; kw...) = open(io -> write_siesta(io, c; kw...), path, "w")

function write_siesta(io::IO, c::Crystal; kspacing::Real=0.25)
    spc = species_list(c)
    title = isempty(c.title) ? formula(c) : c.title
    println(io, "# ", title, " — written by AtomBond")
    println(io, "SystemName ", title)
    println(io, "SystemLabel ", formula(c))
    println(io, "NumberOfAtoms ", natoms(c), "\nNumberOfSpecies ", length(spc))
    println(io, "%block ChemicalSpeciesLabel")
    for (k, s) in enumerate(spc)
        @printf(io, "  %d %d %s\n", k, element(s).Z, s)
    end
    println(io, "%endblock ChemicalSpeciesLabel")
    println(io, "\nLatticeConstant 1.0 Ang\n%block LatticeVectors")
    for v in cell_vectors(c.lattice)
        @printf(io, "  %16.10f %16.10f %16.10f\n", v...)
    end
    println(io, "%endblock LatticeVectors")
    println(io, "\nAtomicCoordinatesFormat Fractional\n%block AtomicCoordinatesAndAtomicSpecies")
    for (f, s) in zip(c.frac, c.species)
        @printf(io, "  %16.10f %16.10f %16.10f %d\n", f..., findfirst(==(s), spc))
    end
    println(io, "%endblock AtomicCoordinatesAndAtomicSpecies")
    n = kpoint_grid(c; kspacing)
    println(io, "\n# k-points: mesh for a spacing of $(kspacing) 1/Å\n%block kgrid_Monkhorst_Pack")
    for k in 1:3
        println(io, "  ", join((j == k ? n[k] : 0 for j in 1:3), " "), " 0.0")
    end
    println(io, "%endblock kgrid_Monkhorst_Pack")
    println(io, "\nPAO.BasisSize DZP\nXC.functional GGA\nXC.authors PBE\nMeshCutoff 300 Ry\nDM.Tolerance 1.0d-5")
    if any(!isnan, c.magmoms)
        println(io, "\nSpin polarized\n%block DM.InitSpin")
        for (i, m) in enumerate(c.magmoms)
            isnan(m) || m == 0 || @printf(io, "  %d %.4f\n", i, m)
        end
        println(io, "%endblock DM.InitSpin")
    end
    has_constraints(c) || return
    cc = c.cell_constraints
    println(io, "\n# Relaxation: ", describe_cell_constraints(cc))
    println(io, "MD.TypeOfRun CG\nMD.NumCGsteps 200\nMD.MaxForceTol 0.02 eV/Ang")
    if cc != (0, 0, 0, 0, 0, 0)
        println(io, "MD.VariableCell T\nMD.MaxStressTol 0.1 GPa")
        cc == FREE_CELL || println(io, "# Note: SIESTA relaxes the whole cell; these cell constraints are not applied")
    end
    full = [i for i in eachindex(c.fixed) if all(c.fixed[i])]
    if !isempty(full)
        println(io, "%block Geometry.Constraints\n  atom ", join(full, " "), "\n%endblock Geometry.Constraints")
    end
    any(f -> any(f) && !all(f), c.fixed) &&
        println(io, "# Note: partly fixed atoms are not written (only fully fixed atoms are constrained)")
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
    ("Quantum ESPRESSO (.in)", :qe),
    ("ABINIT (.abi)", :abinit),
    ("SIESTA (.fdf)", :siesta),
    ("VASP (POSCAR + KPOINTS)", :vasp),
    ("CASTEP (.cell)", :castep),
    ("CIF (P1)", :cif),
    ("Extended XYZ", :xyz),
]

"""
    export_inputs(c, format, dir, name; kspacing=0.25, ecutwfc=40, pseudo="{El}.UPF") -> paths

Write `c` for `format` ∈ (:qe, :abinit, :siesta, :vasp, :castep, :cif, :xyz)
into `dir`. VASP files go into a sub-folder `name`; the others are
`dir/name.<ext>`. `pseudo` is the pseudopotential file pattern for Quantum
ESPRESSO and ABINIT (`{El}` becomes the element; ABINIT gets `.psp8` when
the pattern is the default `.UPF` one).
"""
function export_inputs(c::Crystal, format::Symbol, dir::AbstractString, name::AbstractString;
                       kspacing::Real=0.25, ecutwfc::Real=40, pseudo::AbstractString="{El}.UPF")
    isempty(strip(name)) && throw(ArgumentError("empty file name"))
    mkpath(dir)
    format === :vasp && return write_vasp_inputs(joinpath(dir, name), c; kspacing)
    ext = Dict(:qe => ".in", :abinit => ".abi", :siesta => ".fdf", :castep => ".cell", :cif => ".cif", :xyz => ".xyz")
    haskey(ext, format) || throw(ArgumentError("unknown export format $format"))
    path = joinpath(dir, name * ext[format])
    format === :qe ? write_qe(path, c; kspacing, ecutwfc, pseudo) :
    format === :abinit ? write_abinit(path, c; kspacing, ecutwfc,
                                      pseudo=pseudo == "{El}.UPF" ? "{El}.psp8" : pseudo) :
    format === :siesta ? write_siesta(path, c; kspacing) :
    format === :castep ? write_castep(path, c; kspacing) :
    format === :cif ? write_cif(path, c) : write_xyz(path, c)
    return [path]
end
