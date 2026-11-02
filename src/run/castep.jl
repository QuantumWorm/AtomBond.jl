# CASTEP as a DFT engine: `.cell` from write_castep, a `.param` for the task,
# results from the `.castep` output and the `.bands` file, relaxation steps
# from the `.geom` file. CASTEP is licensed software: AtomBond only finds and
# runs the user's own copy.

"""
CASTEP. Found as `castep.serial` or `castep` on the PATH (e.g. in `~/bin`),
or by the `command` in the settings file (e.g. `"mpirun -np 4 castep.mpi"`);
the seed name is added. Pseudopotentials: empty (or `on the fly`) lets
CASTEP generate them; a folder is searched for a `.usp`, `.recpot` or `.upf`
file per element (name starting with the element symbol), which is copied
next to the input.
"""
struct CASTEPEngine <: DFTEngine end
engine_id(::CASTEPEngine) = :castep
engine_name(::CASTEPEngine) = "CASTEP"
engine_executables(::CASTEPEngine) = ["castep.serial", "castep"]
engine_install_hint(::CASTEPEngine) = "licensed: put castep in ~/bin or set its command"
engine_pseudo_default(e::CASTEPEngine) = something(engine_setting(e, "pseudo_dir"), "on the fly")
engine_command(::CASTEPEngine, command::Vector{String}, dir, name, input) = `$command $name`
function engine_finished(::CASTEPEngine, job)
    isfile(job.results) || return false
    return any(l -> occursin("Total time", l), readlines(job.results))
end

"""
    write_castep_param(io; task="SinglePoint", ecutwfc=40, stress=false, spin=false)

A CASTEP `.param` file for a run AtomBond reads back: PBE, `cut_off_energy`
from `ecutwfc` (Ry, as the other codes; written in eV), Gaussian smearing of
1 mHa (as DFTK), extra empty bands (for the gap), the stress when the cell
relaxes, and geometry tolerances of 5 meV/Å and 0.05 GPa.
"""
function write_castep_param(io::IO; task::AbstractString="SinglePoint", ecutwfc::Real=40, stress::Bool=false,
                            spin::Bool=false)
    println(io, "! written by AtomBond")
    println(io, "task                : ", task)
    println(io, "xc_functional       : PBE")
    @printf(io, "cut_off_energy      : %.2f eV   ! %.2f Ry\n", ecutwfc * RY_EV, ecutwfc)
    println(io, "elec_energy_tol     : 1.0e-7 eV")
    println(io, "fix_occupancy       : false")
    println(io, "smearing_scheme     : gaussian")
    println(io, "smearing_width      : 0.0272 eV")
    println(io, "nextra_bands        : 8")
    spin && println(io, "spin_polarized      : true")
    if task == "GeometryOptimization"
        println(io, "geom_force_tol      : 0.005 eV/ang")
        println(io, "geom_max_iter       : 100")
        stress && println(io, "geom_stress_tol     : 0.05 GPa\ncalculate_stress    : true\nfinite_basis_corr   : 2")
    end
    println(io, "write_checkpoint    : none")
end

function write_engine_inputs(::CASTEPEngine, dir, name, c::Crystal; relax::Symbol, kspacing, ecutwfc, pseudo)
    cc = relax === :atoms ? with_constraints(c; cell=(0, 0, 0, 0, 0, 0)) :
         relax === :cell && c.cell_constraints == (0, 0, 0, 0, 0, 0) ? with_constraints(c; cell=FREE_CELL) : c
    input = joinpath(dir, name * ".cell")
    open(input, "w") do io
        write_castep(io, cc; kspacing)
        p = strip(pseudo)
        if !(isempty(p) || occursin(r"^on.the.fly$"i, p))
            pdir = abspath(expanduser(p))
            isdir(pdir) || throw(ArgumentError("no pseudopotential folder $pdir (or use 'on the fly')"))
            files = filter(f -> occursin(r"\.(usp|recpot|upf)$"i, f), readdir(pdir))
            println(io, "\n%BLOCK SPECIES_POT")
            for el in species_list(c)
                m = filter(f -> occursin(Regex("^" * string(el) * "[._-]", "i"), f), files)
                isempty(m) && throw(ArgumentError("no pseudopotential for $el in $pdir"))
                f = sort(m; by=f -> (!startswith(f, string(el)), length(f)))[1]
                cp(joinpath(pdir, f), joinpath(dir, f); force=true)
                println(io, "  ", rpad(string(el), 3), " ", f)
            end
            println(io, "%ENDBLOCK SPECIES_POT")
        end
    end
    open(joinpath(dir, name * ".param"), "w") do io
        write_castep_param(io; task=relax === :none ? "SinglePoint" : "GeometryOptimization", ecutwfc,
                           stress=relax === :cell, spin=any(!isnan, c.magmoms))
    end
    return input, joinpath(dir, name * ".castep")
end

# --- reading ------------------------------------------------------------------------

"""
    read_castep_bands(path) -> (; eigenvalues, kweights, fermi, n_spin)

A CASTEP `.bands` file: the levels (eV) per k-point, spin up for every
k-point first (as DFTK), the k-point weights (sum 1) and the Fermi level (eV,
the first spin's when there are two).
"""
function read_castep_bands(path::AbstractString)
    lines = readlines(path)
    val(key) = split(lines[findfirst(l -> startswith(l, key), lines)])
    nk = parse(Int, val("Number of k-points")[end])
    nspin = parse(Int, val("Number of spin components")[end])
    fl = findfirst(l -> startswith(l, "Fermi energ"), lines)
    fermi = HARTREE_EV * parse(Float64, split(lines[fl])[end-nspin+1])
    neig = [parse(Int, w) for w in val("Number of eigenvalues")[end-nspin+1:end]]
    eigs = [[Vector{Float64}() for _ in 1:nk] for _ in 1:nspin]
    w = zeros(nk)
    k, σ = 0, 0
    for l in lines
        t = split(l)
        isempty(t) && continue
        if t[1] == "K-point"
            k = parse(Int, t[2]); w[k] = parse(Float64, t[6]); σ = 0
        elseif t[1] == "Spin" && k > 0
            σ = parse(Int, t[3])
        elseif k > 0 && σ > 0 && length(t) == 1 && (x = tryparse(Float64, t[1])) !== nothing
            length(eigs[σ][k]) < neig[σ] && push!(eigs[σ][k], HARTREE_EV * x)
        end
    end
    return (; eigenvalues=[e for s in eigs for e in s], kweights=repeat(w ./ sum(w), nspin), fermi, n_spin=nspin)
end

"""Rows of the last starred block (forces, stress) after `header` in `lines`: the numbers of each row."""
function castep_block(lines, header::Regex)
    k = findlast(l -> occursin(header, l), lines)
    k === nothing && return nothing
    rows = Vector{Vector{String}}()
    for l in lines[k+1:end]
        startswith(strip(l), "****") && break
        t = split(strip(strip(l), '*'))
        length(t) >= 4 && all(x -> tryparse(Float64, x) !== nothing, t[end-2:end]) && push!(rows, t)
    end
    return rows
end

"""
    read_castep_results(path; title="", clamped=nothing) -> DFTResults

The results of a CASTEP run from its `.castep` output: the final structure
(as [`read_castep_output`](@ref)), the final free energy, forces and stress,
and, from the `.bands` file next to it, the levels and the Fermi level; for a
geometry optimization every step from the `.geom` file.
"""
function read_castep_results(path::AbstractString; title::AbstractString="", clamped=nothing)
    lines = readlines(path)
    c = read_castep_output(path; charges=:none)
    c = Crystal(c.lattice, c.species, c.frac; pbc=c.pbc, title)
    efree = findlast(l -> occursin(r"Final free energy \(E-TS\)", l), lines)
    efin = findlast(l -> occursin(r"Final energy(, E)?\s*=", l), lines)
    k = something(efree, efin, 0)
    k == 0 && throw(ArgumentError("no final energy in $path (CASTEP did not finish?)"))
    energy = parse(Float64, match(r"=\s*([-\d.Ee+]+)\s*eV", lines[k])[1])
    frows = castep_block(lines, r"\*+ .*Forces \*+")
    forces = frows === nothing ? Vec3[] : [Vec3(parse.(Float64, r[end-2:end])) for r in frows]
    geom = any(l -> occursin(r"task\s*:\s*GeometryOptimi[sz]ation"i, l), lines) ||
           any(l -> occursin("BFGS: finished iteration", l) || occursin("LBFGS: finished iteration", l), lines)
    conv = geom ? any(l -> occursin(r"Geometry optimi[sz]ation completed successfully", l), lines) :
           !any(l -> occursin(r"did not converge|SCF cycles failed"i, l), lines)
    its = [parse(Int, m[1]) for m in (match(r"BFGS: finished iteration\s+(\d+)", l) for l in lines) if m !== nothing]
    steps = isempty(its) ? 0 : maximum(its)              # iteration 0 is the starting structure
    scf_start = findlast(l -> occursin(r"^SCF loop", l), lines)
    nscf = scf_start === nothing ? 0 :
           count(l -> occursin(r"^\s+\d+\s+-?[\d.]+E[-+]\d+.*<-- SCF", l), lines[scf_start:end])
    vline = findfirst(l -> occursin(r"CASTEP version", l), lines)
    version = vline === nothing ? "" : something(match(r"version\s+([\d.]+)", lines[vline]), (; captures=[""])).captures[1]
    seed = path[1:end-length(".castep")]
    bands = isfile(seed * ".bands") ? read_castep_bands(seed * ".bands") : nothing
    eigs = bands === nothing ? Vector{Float64}[] : bands.eigenvalues
    fermi = bands === nothing ? nothing : bands.fermi
    nspin = bands === nothing ? 1 : bands.n_spin
    occs = [Float64[e <= fermi + 1e-6 ? 2 / nspin : 0.0 for e in es] for es in eigs]   # .bands has no occupations
    history = geom && isfile(seed * ".geom") ? read_castep_history(seed * ".geom", c; clamped) : nothing
    return DFTResults(; calculation=geom ? "relax" : "scf", converged=conv, iterations=geom ? steps : nscf, energy,
                      structure=c, fermi, forces, eigenvalues=eigs, occupations=occs,
                      kweights=bands === nothing ? Float64[] : bands.kweights, n_spin=nspin, version, history,
                      program="CASTEP")
end

"""
    read_castep_history(geom, like; clamped=nothing) -> RelaxHistory or nothing

The steps of a CASTEP geometry optimization from its `.geom` file (atomic
units): energy, cell, positions, forces and, when computed, the stress
(CASTEP's sign: positive when tensile, as AtomBond's). Species and title
come from `like`.
"""
function read_castep_history(path::AbstractString, like::Crystal; clamped=nothing)
    frames = NamedTuple[]
    cur = nothing
    function finish!()
        cur === nothing || isempty(cur.R) || push!(frames, cur)
        cur = nothing
    end
    for raw in eachline(path)
        parts = split(raw, "<--")
        length(parts) < 2 && continue
        tag, w = strip(parts[2]), split(parts[1])
        if tag == "E"
            finish!()
            cur = (; E=parse(Float64, w[1]), h=Vec3[], S=Vector{Float64}[], R=Vec3[], F=Vec3[], sp=Symbol[])
        elseif cur !== nothing
            tag == "h" && push!(cur.h, BOHR * Vec3(parse.(Float64, w[1:3])))
            tag == "S" && push!(cur.S, parse.(Float64, w[1:3]))
            tag == "R" && (push!(cur.R, BOHR * Vec3(parse.(Float64, w[3:5]))); push!(cur.sp, castep_species(w[1])))
            tag == "F" && push!(cur.F, HA_PER_BOHR_EV_PER_A * Vec3(parse.(Float64, w[3:5])))
        end
    end
    finish!()
    frames = filter(f -> length(f.h) == 3 && length(f.F) == length(f.R), frames)
    isempty(frames) && return nothing
    n = length(frames[1].R)
    structures = [crystal_from_cartesian(Lattice(f.h...), f.sp, f.R; pbc=like.pbc,
                                         title=like.title * " (CASTEP step $(k - 1))") for (k, f) in enumerate(frames)]
    cl = clamped === nothing || length(clamped) != n ? fill(false, n) : collect(Bool, clamped)
    free = findall(!, cl)
    maxf = [isempty(free) ? 0.0 : maximum(norm, f.F[free]) for f in frames]
    to_gpa = HARTREE_EV / BOHR^3 * EV_PER_A3_IN_GPA
    stresses = all(f -> length(f.S) == 3, frames) && any(f -> any(r -> any(!iszero, r), f.S), frames) ?
               [to_gpa .* permutedims(reduce(hcat, f.S)) for f in frames] : Matrix{Float64}[]
    return RelaxHistory([HARTREE_EV * f.E for f in frames], maxf, structures, [f.F for f in frames], cl, stresses)
end

function read_engine_results(::CASTEPEngine, job)
    like, clamped = grouped_like(job.crystal)
    return read_castep_results(job.results; title=like.title, clamped)
end

function engine_history(::CASTEPEngine, job)
    geom = job.results[1:end-length(".castep")] * ".geom"
    job.relax && isfile(geom) || return nothing
    like, clamped = grouped_like(job.crystal)
    return try
        read_castep_history(geom, like; clamped)
    catch
        nothing
    end
end

register_engine!(CASTEPEngine())
