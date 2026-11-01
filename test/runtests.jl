using Test
ENV["ATOMBOND_HOME"] = mktempdir()      # recent-files list of the tests, not the user's
using AtomBond
using LinearAlgebra
using AtomBond: RGBf

const STRUCTS = joinpath(@__DIR__, "..", "examples", "structures")

@testset "AtomBond" begin
    @testset "elements" begin
        @test element(:Si).Z == 14
        @test element("SI").symbol == :Si
        @test element("Si1").symbol == :Si
        @test element("O2-").symbol == :O
        @test element(8).symbol == :O
        @test element("Xx").Z == 0          # unknown → placeholder, no throw
        # The table: Z 1–103 with sources' values
        @test length(AtomBond.ELEMENTS) == 103 && all(z -> element(z).Z == z, 1:103)
        @test element(:Ce).name == "Cerium" && element(58).symbol == :Ce
        @test element(:O).color == AtomBond.hexcolor("FF0D0D") && element(:N).color == AtomBond.hexcolor("3050F8")
        @test element(:H).covalent_radius == 0.31 < element(:Si).covalent_radius == 1.11
        @test atomic_mass(:Si) == 28.085 && atomic_mass(:Be) == 9.0121831 && atomic_mass(:Tc) == 98
        @test_throws ArgumentError atomic_mass(:X)
        @test element_property(:Si, :vdw_radius_alvarez) == 2.19 && element_property(:Si, :vdw_radius_bondi) == 2.10
        @test element_property(:Lr, :covalent_radius) === nothing && element_property(:Pu, :mass_note) == "mass number"
        @test vdw_radius(:Si) == 2.10 && vdw_radius(:Fe) == 2.44
        # Overrides: a user table, and VESTA's elements.ini
        saved = deepcopy(AtomBond.ELEMENT_TABLE)
        try
            dir = mktempdir()
            write(joinpath(dir, "elements.csv"), "# mine\nsymbol,color,vdw_radius,electronegativity\nO,00FF00,1.4,3.44\n")
            override_elements!(joinpath(dir, "elements.csv"))
            @test element(:O).color == AtomBond.hexcolor("00FF00") && vdw_radius(:O) == 1.4
            @test element_property(:O, :electronegativity) == 3.44 && atomic_mass(:O) == 15.999
            write(joinpath(dir, "elements.ini"), "  1  H  0.46  1.20  0.200    1.00000    0.80000    0.80000\n" *
                                                 "  1  D  0.46  1.20  0.200    0.80000    0.80000    1.00000\n" *
                                                 " 14 Si  1.18  2.10  0.260    0.10596    0.23226    0.98096\n")
            @test use_vesta_colors!(joinpath(dir, "elements.ini")) == 2
            @test element(:Si).color == AtomBond.hexcolor("1B3BFA") && element(:H).color == AtomBond.hexcolor("FFCCCC")
            withenv("ATOMBOND_HOME" => dir) do
                write(joinpath(dir, "elements.csv"), "symbol,color\nSi,123456\n")
                @test AtomBond.load_user_elements!() == joinpath(dir, "elements.csv")
            end
            @test element(:Si).color == AtomBond.hexcolor("123456")
        finally
            empty!(AtomBond.ELEMENT_TABLE); merge!(AtomBond.ELEMENT_TABLE, saved); AtomBond.refresh_elements!()
        end
        @test element(:O).color == AtomBond.hexcolor("FF0D0D")
    end

    @testset "lattice" begin
        L = Lattice(3.0, 4.0, 5.0, 80, 95, 110)
        p = lattice_parameters(L)
        @test all(isapprox.(p, (3.0, 4.0, 5.0, 80.0, 95.0, 110.0); atol=1e-10))
        f = [0.1, 0.7, 0.35]
        @test cart_to_frac(L, frac_to_cart(L, f)) ≈ f
        @test cell_volume(cubic(2.0)) ≈ 8.0
        @test_throws ArgumentError Lattice([1 0 0; 0 1 0; 0 0 0])
    end

    @testset "crystal" begin
        c = Crystal(cubic(5.0), ["Ga", "As"], [[0, 0, 0], [0.25, 0.25, 1.25]])
        @test natoms(c) == 2
        @test formula(c) == "Ga1As1"
        @test wrap(c).frac[2] ≈ [0.25, 0.25, 0.25]
        @test_throws ArgumentError Crystal(cubic(5.0), [:Si], [[0, 0, 0], [0.5, 0.5, 0.5]])
    end

    @testset "POSCAR" begin
        si = read_poscar(joinpath(STRUCTS, "Si_conv.vasp"))
        @test natoms(si) == 8
        @test cell_volume(si.lattice) ≈ 5.431^3
        # Nearest-neighbour distance in diamond is a√3/4.
        r = cartesian_positions(si)
        @test norm(r[5] - r[1]) ≈ 5.431 * sqrt(3) / 4

        gan = read_poscar(joinpath(STRUCTS, "GaN_wurtzite.vasp"))
        a, b, c, α, β, γ = lattice_parameters(gan.lattice)
        @test γ ≈ 120 atol = 1e-3
        @test species_list(gan) == [:Ga, :N]

        # Round trip, both Direct and Cartesian.
        for cart in (false, true)
            io = IOBuffer()
            write_poscar(io, gan; cartesian=cart)
            back = read_poscar(IOBuffer(String(take!(io))))
            @test back.species == gan.species
            @test all(isapprox.(back.frac, gan.frac; atol=1e-10))
            @test back.lattice.matrix ≈ gan.lattice.matrix
        end

        # Negative scale = target volume.
        txt = replace(read(joinpath(STRUCTS, "Si_conv.vasp"), String), "\n5.431\n" => "\n-27.0\n")
        @test cell_volume(read_poscar(IOBuffer(txt)).lattice) ≈ 27.0
    end

    @testset "XYZ" begin
        w = read_xyz(joinpath(STRUCTS, "water.xyz"))
        @test w.species == [:O, :H, :H]
        @test w.pbc == (false, false, false)
        r = cartesian_positions(w)
        @test norm(r[2] - r[1]) ≈ 0.9572 atol = 1e-3

        # extxyz round trip keeps the cell.
        si = load_structure(joinpath(STRUCTS, "Si_conv.vasp"))
        io = IOBuffer()
        write_xyz(io, si)
        back = read_xyz(IOBuffer(String(take!(io))))
        @test back.lattice.matrix ≈ si.lattice.matrix
        @test all(isapprox.(back.frac, si.frac; atol=1e-9))
        @test back.pbc == (true, true, true)
    end

    @testset "bonds" begin
        load(f) = load_structure(joinpath(STRUCTS, f))
        # Expected bond counts: every atom 4-fold (diamond, zincblende,
        # wurtzite) or 6-fold (rocksalt) coordinated → N*z/2 bonds.
        for (f, nb, z) in [("Si_conv.vasp", 16, 4), ("Si_prim.vasp", 4, 4),
                           ("GaAs_conv.vasp", 16, 4), ("GaN_wurtzite.vasp", 8, 4),
                           ("NaCl_conv.vasp", 24, 6), ("water.xyz", 2, nothing)]
            c = load(f)
            b = find_bonds(c)
            @test length(b) == nb
            z === nothing || @test all(==(z), coordination(c, b))
        end
        si = load("Si_conv.vasp")
        @test all(b -> b.length ≈ 5.431 * sqrt(3) / 4, find_bonds(si))
        # A cell much smaller than the cutoff still finds every image.
        @test length(find_bonds(Crystal(cubic(2.3), [:Si], [[0, 0, 0]]))) == 3
        # Bonds don't depend on where atoms sit relative to the cell.
        shifted = Crystal(si.lattice, si.species, [f .+ [0.37, -1.2, 2.05] for f in si.frac])
        @test length(find_bonds(shifted)) == 16
    end

    @testset "display model" begin
        si = load_structure(joinpath(STRUCTS, "Si_conv.vasp"))
        dm = AtomBond.build_display(si, ViewSettings())
        # 8 corners + 6 face centres + 4 interior atoms.
        @test length(dm.atoms) == 18
        # Every bond of the 4 interior atoms ends on a drawn atom, so those
        # 16 bonds are whole (2 halves each); the remaining halves are stubs.
        full = 2 * 16
        @test length(dm.halfbonds) >= full
        nostub = AtomBond.build_display(si, ViewSettings(bond_stubs=false))
        @test length(nostub.halfbonds) == full
        @test length(AtomBond.build_display(si, ViewSettings(show_boundary_atoms=false)).atoms) == 8
        @test isempty(AtomBond.build_display(si, ViewSettings(show_bonds=false)).halfbonds)

        w = AtomBond.build_display(load_structure(joinpath(STRUCTS, "water.xyz")), ViewSettings())
        @test length(w.atoms) == 3          # no images for a molecule
        @test length(w.halfbonds) == 4
    end

    @testset "display styles" begin
        si = load_structure(joinpath(STRUCTS, "Si_conv.vasp"))
        r(style) = AtomBond.atom_radius(:Si, ViewSettings(; style))
        @test r(:cpk) ≈ vdw_radius(:Si) ≈ 2.10
        @test r(:stick) ≈ ViewSettings().bond_radius
        @test r(:ball_and_stick) ≈ 0.4 * 1.11
        @test vdw_radius(:Mo) ≈ 2.45                                 # Alvarez where Bondi has none
        @test vdw_radius(:Fm) ≈ element(:Fm).covalent_radius + 0.8   # neither: fallback
        # CPK has no bonds; the other styles do.
        @test isempty(AtomBond.build_display(si, ViewSettings(style=:cpk)).bonds)
        for s in (:ball_and_stick, :stick, :line)
            @test length(AtomBond.build_display(si, ViewSettings(style=s)).bonds) == 16
        end
    end

    @testset "view presets" begin
        gan = load_structure(joinpath(STRUCTS, "GaN_wurtzite.vasp"))
        L = gan.lattice
        a, b, c = cell_vectors(L)
        @test view_direction(L, :c)[1] ≈ c
        @test view_direction(L, [1, 1, 0])[1] ≈ a + b
        # c* is perpendicular to a and b; a* to b and c.
        cs = view_direction(L, :cstar)[1]
        @test abs(dot(cs, a)) < 1e-10 && abs(dot(cs, b)) < 1e-10
        as = view_direction(L, :astar)[1]
        @test abs(dot(as, b)) < 1e-10 && abs(dot(as, c)) < 1e-10
        @test_throws ArgumentError view_direction(L, :d)
    end

    @testset "CIF" begin
        W, w = AtomBond.parse_symop("-y+1/2, x-y, 'z+2/3'")
        @test W == [0 -1 0; 1 -1 0; 0 0 1]
        @test w ≈ [0.5, 0, 2 / 3]
        @test AtomBond.parse_symop("1/2+X,-Y,0.5-z")[2] ≈ [0.5, 0, 0.5]
        @test AtomBond.cif_number("5.4310(12)") == 5.431
        @test AtomBond.hall_number_from_hm("F d -3 m :2") == 526
        @test AtomBond.hall_number_from_hm("F d -3 m Z") == 526      # ICSD origin-2 suffix
        @test AtomBond.hall_number_from_hm("R -3 m") == 458        # hexagonal axes by default
        @test_throws ErrorException AtomBond.hall_number_from_hm("Q 9 z")

        # (file, atoms, space group, coordination per species)
        for (f, n, sg, z) in [("TiO2_rutile.cif", 6, 136, Dict(:Ti => 6, :O => 3)),     # explicit ops
                              ("ZnO_wurtzite.cif", 4, 186, Dict(:Zn => 4, :O => 4)),    # explicit, unquoted
                              ("Si_Fd-3m.cif", 8, 227, Dict(:Si => 4)),                 # H-M symbol only
                              ("GaN_P63mc.cif", 4, 186, Dict(:Ga => 4, :N => 4))]       # hexagonal, symbol only
            c = load_structure(joinpath(STRUCTS, f))
            @test natoms(c) == n
            @test spacegroup(c)[1] == sg
            zc = coordination(c, find_bonds(c))
            @test all(zc[i] == z[c.species[i]] for i in 1:n)
        end
        @test load_structure(joinpath(STRUCTS, "TiO2_rutile.cif")).title == "Rutile"
        # Same structure from POSCAR and from CIF.
        @test cell_volume(load_structure(joinpath(STRUCTS, "Si_Fd-3m.cif")).lattice) ≈ 5.431^3
    end

    @testset "bonding rules" begin
        @test ismetal(:Ti) && !ismetal(:O) && !ismetal(:Si)
        rut = load_structure(joinpath(STRUCTS, "TiO2_rutile.cif"))
        # Ti–Ti along c (2.96 Å) is under the covalent cutoff; only kept on request.
        @test length(find_bonds(rut; metal_metal=true, bridge_angle=nothing)) > length(find_bonds(rut))
        # A pure metal keeps its metal–metal bonds automatically.
        cu = supercell(Crystal(cubic(3.615), [:Cu, :Cu, :Cu, :Cu],
                               [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]), 1, 1, 1)
        @test all(==(12), coordination(cu, find_bonds(cu)))
    end

    @testset "symmetry and supercells" begin
        si = load_structure(joinpath(STRUCTS, "Si_conv.vasp"))
        @test spacegroup(si) == (227, "Fd-3m")
        p = primitive_cell(si)
        @test natoms(p) == 2
        @test cell_volume(p.lattice) ≈ cell_volume(si.lattice) / 4
        @test natoms(conventional_cell(p)) == 8
        @test spacegroup(load_structure(joinpath(STRUCTS, "water.xyz"))) === nothing

        s = supercell(si, 2, 1, 3)
        @test natoms(s) == 48
        @test cell_volume(s.lattice) ≈ 6 * cell_volume(si.lattice)
        @test spacegroup(s)[1] == 227
        @test length(find_bonds(s)) == 6 * 16
        # General matrix: √2×√2 rotated cell.
        r = supercell(si, [1 -1 0; 1 1 0; 0 0 1])
        @test natoms(r) == 16
        @test all(==(4), coordination(r, find_bonds(r)))
        @test_throws ArgumentError supercell(si, [1 0 0; 0 1 0; 0 0 0])
        # Molecules can't be repeated along their non-periodic box.
        @test_throws ArgumentError supercell(load_structure(joinpath(STRUCTS, "water.xyz")), 2, 1, 1)
    end

    @testset "geometry" begin
        @test distance([0, 0, 0], [3, 4, 0]) ≈ 5
        @test bond_angle([1, 0, 0], [0, 0, 0], [0, 1, 0]) ≈ 90
        @test bond_angle([1, 1, 1], [0, 0, 0], [1, -1, -1]) ≈ acosd(-1 / 3)    # tetrahedral
        # Dihedral sign: looking down q→r, clockwise rotation of p onto s is positive.
        p, q, r = [1.0, 0, 0], [0.0, 0, 0], [0.0, 0, 1]
        @test dihedral(p, q, r, [0.0, 1, 1]) ≈ 90
        @test dihedral(p, q, r, [0.0, -1, 1]) ≈ -90
        @test abs(dihedral(p, q, r, [-1.0, 0, 1])) ≈ 180
        @test measure([p, q]) == (:distance, 1.0)
        @test measure([p]) === nothing

        si = load_structure(joinpath(STRUCTS, "Si_conv.vasp"))
        nn = 5.431 * sqrt(3) / 4
        @test min_image_distance(si, 1, 5) ≈ nn
        # Atom 8 at (¾,¾,¼) is a nearest neighbour of atom 1 only through the boundary.
        @test min_image_distance(si, 1, 8) ≈ nn
        @test min_image_distance(si, 1, 2) ≈ 5.431 / sqrt(2)
    end

    @testset "selection" begin
        v = visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false, panel=false)
        @test v.info[] == AtomBond.SELECTION_HINT
        select!(v, 5)
        @test startswith(v.info[], "Si5") && occursin("4 bonds", v.info[])
        select!(v, 1; add=true)
        @test occursin("distance 2.3517 Å", v.info[])
        select!(v, 2; add=true)                         # Si5 – Si1 – Si2: not bonded at Si1
        @test occursin("angle", v.info[])
        select!(v, 2; add=true)                         # Shift+click again removes it
        @test length(v.selection) == 2
        @test length(v.overlay) == 3                    # halo, dashed line, value
        # An image atom that disappears with the boundary atoms: its own copy is selected.
        select!(v, 1, (1, 1, 1))
        v.settings.show_boundary_atoms = false
        redraw!(v)
        @test v.selection == [(1, AtomBond.Shift(0, 0, 0))]
        v.settings.show_boundary_atoms = true
        select!(v, 3)
        set_crystal!(v, supercell(v.crystal, 2, 1, 1))  # new structure: selection reset
        @test isempty(v.selection) && isempty(v.overlay)
        clear_selection!(v)
        @test v.info[] == AtomBond.SELECTION_HINT
    end

    @testset "editing" begin
        si = load_structure(joinpath(STRUCTS, "Si_conv.vasp"))
        d = delete_atoms(si, [1, 5])
        @test natoms(d) == 6 && natoms(si) == 8            # input untouched
        @test d.frac[1] == si.frac[2]
        @test_throws ArgumentError delete_atoms(si, [9])
        g = substitute(si, [5], "ge")
        @test g.species[5] == :Ge && formula(g) == "Si7Ge1"
        @test_throws ArgumentError substitute(si, [1], "Qx")
        # Interstitial at the tetrahedral site (½,½,½) is 2.35 Å from its neighbours.
        t = add_atom(si, :Si, [0.5, 0.5, 0.5])
        @test natoms(t) == 9
        @test_throws ArgumentError add_atom(si, :Si, [0.01, 0, 0])           # on top of Si1
        @test_throws ArgumentError add_atom(si, :Si, [0.999, 0, 0])          # same, across the boundary
        @test natoms(add_atom(si, :O, [2.7, 0, 0]; cartesian=true)) == 9
        m = move_atoms(si, [5], [0.1, 0, 0])
        @test cartesian_positions(m)[5] ≈ cartesian_positions(si)[5] + [0.1, 0, 0]
    end

    @testset "undo / redo" begin
        v = visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false, panel=false)
        select!(v, 5); select!(v, 6; add=true)
        AtomBond.substitute_selected!(v, :Ge)
        @test formula(v.crystal) == "Si6Ge2"
        @test length(v.selection) == 2                     # kept: numbering unchanged
        @test length(v.legend.entrygroups[][1][2]) == 2    # legend now lists Si and Ge
        AtomBond.delete_selected!(v)
        @test natoms(v.crystal) == 6 && isempty(v.selection)
        @test undo!(v) && formula(v.crystal) == "Si6Ge2"
        @test undo!(v) && formula(v.crystal) == "Si8"
        @test !undo!(v) && v.status[] == "Nothing to undo"
        @test redo!(v) && formula(v.crystal) == "Si6Ge2"
        # A new edit clears the redo stack.
        AtomBond.edit!(c -> supercell(c, 2, 1, 1), v, "make supercell"; reframe=true)
        @test natoms(v.crystal) == 16 && isempty(v.future)
        # A failing edit changes nothing and reports why.
        @test !AtomBond.edit!(c -> add_atom(c, :Si, [0, 0, 0]), v, "add Si")
        @test natoms(v.crystal) == 16 && startswith(v.status[], "Could not add Si")
        clear_selection!(v)
        AtomBond.delete_selected!(v)
        @test occursin("Select atoms first", v.status[])

        # Keyboard: Delete, Ctrl+Z, Ctrl+Y (key events as the window would send them).
        M = AtomBond.Makie
        ev = M.events(v.lscene.scene)
        key(k, a=M.Keyboard.press) = (ev.keyboardbutton[] = M.KeyEvent(k, a))
        select!(v, 1)
        key(M.Keyboard.delete)
        @test natoms(v.crystal) == 15
        key(M.Keyboard.left_control); key(M.Keyboard.z)
        @test natoms(v.crystal) == 16
        key(M.Keyboard.y)
        @test natoms(v.crystal) == 15
        key(M.Keyboard.left_control, M.Keyboard.release)
    end

    @testset "surfaces" begin
        si = load_structure(joinpath(STRUCTS, "Si_conv.vasp"))
        sp = primitive_cell(si)
        a0 = 5.431
        # (crystal, plane, expected interplanar spacing of the cell)
        for (c, hkl, d) in [(si, (0, 0, 1), a0), (si, (1, 1, 0), a0 / sqrt(2)), (si, (1, 1, 1), a0 / sqrt(3)),
                            (sp, (1, 1, 1), a0 / sqrt(3)), (si, (2, 1, 3), a0 / sqrt(14)),
                            (si, (1, -1, 2), a0 / sqrt(6)), (si, (0, 0, -1), a0), (si, (2, 2, 2), a0 / sqrt(3))]
            P = AtomBond.surface_basis(c.lattice, hkl)
            @test round(Int, det(P)) == 1
            @test dot(collect(hkl) .÷ gcd(hkl...), P[:, 3]) == 1
            s = surface_cell(c, hkl)
            a, b, _ = cell_vectors(s.lattice)
            @test natoms(s) == natoms(c)
            @test abs(a[3]) < 1e-9 && abs(b[3]) < 1e-9 && abs(a[2]) < 1e-9   # a ∥ x, b in xy
            @test 90 - 1e-6 <= lattice_parameters(s.lattice)[6] <= 120 + 1e-6
            @test AtomBond.layer_spacing(s.lattice) ≈ d
            @test all(==(4), coordination(s, find_bonds(s)))                 # still bulk
        end
        @test_throws ArgumentError surface_cell(si, (0, 0, 0))

        # Si(001) slab, 2 layers, 10 Å vacuum: 16 atoms, c normal to the surface,
        # the two outermost atoms on each side have 2 dangling bonds.
        sl = slab(si, (0, 0, 1); layers=2, vacuum=10, shift=0)
        @test natoms(sl) == 16
        c3 = cell_vectors(sl.lattice)[3]
        @test c3[1] ≈ 0 atol = 1e-9
        @test c3[2] ≈ 0 atol = 1e-9
        @test c3[3] ≈ 2a0 - a0 / 4 + 10
        z = coordination(sl, find_bonds(sl))
        @test count(==(2), z) == 4 && count(==(4), z) == 12
        zf = [f[3] for f in sl.frac]
        @test minimum(zf) * c3[3] ≈ 5 && (1 - maximum(zf)) * c3[3] ≈ 5   # centred
        # Termination shift by a quarter layer moves the cut to the next atomic plane.
        @test slab(si, (0, 0, 1); layers=2, vacuum=10, shift=0.25).frac != sl.frac
        @test natoms(slab(sp, (1, 1, 1); layers=3, vacuum=0)) == 6
        # Automatic termination: Si(111) cut between bilayers → every surface atom
        # keeps 3 of its 4 bonds; Si(001) → 2 dangling bonds per surface atom.
        s111 = slab(si, (1, 1, 1); layers=2, vacuum=10)
        z = coordination(s111, find_bonds(s111))
        @test sort(unique(z)) == [3, 4]
        s001 = slab(si, (0, 0, 1); layers=2, vacuum=10)
        @test sort(unique(coordination(s001, find_bonds(s001)))) == [2, 4]
        @test_throws ArgumentError slab(si, (0, 0, 1); layers=0)
    end

    @testset "superlattices" begin
        si = load_structure(joinpath(STRUCTS, "Si_conv.vasp"))
        ge = Crystal(cubic(5.658), fill(:Ge, 8), si.frac)
        s001 = surface_cell(si, (0, 0, 1))
        # Stacking a cell on itself = 1×1×2 supercell.
        st = stack_cells(s001, s001)
        @test natoms(st) == 16
        @test cell_volume(st.lattice) ≈ 2 * cell_volume(si.lattice)
        @test all(==(4), coordination(st, find_bonds(st)))
        # Si/Ge needs strain (4 % mismatch).
        @test_throws ArgumentError superlattice(si, ge, (0, 0, 1); layers=(2, 2))
        sl = superlattice(si, ge, (0, 0, 1); layers=(2, 3), strain=true)
        @test formula(sl) == "Si16Ge24"
        a, b, c = lattice_parameters(sl.lattice)
        @test a ≈ 5.431 && b ≈ 5.431 && c ≈ 2 * 5.431 + 3 * 5.658
        @test all(==(4), coordination(sl, find_bonds(sl)))                # continuous crystal
        sl111 = superlattice(si, si, (1, 1, 1); layers=(1, 1))
        @test natoms(sl111) == 16
    end

    @testset "surface dialog" begin
        v = visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false, panel=false)
        d = AtomBond.surface_dialog(v)
        d.hkl.stored_string[] = "1 1 0"
        d.layers.stored_string[] = "2"
        d.vacuum.stored_string[] = "12"
        d.build_slab.clicks[] += 1
        @test natoms(v.crystal) == 16 && occursin("(1 1 0) slab", v.status[])
        undo!(v)
        @test natoms(v.crystal) == 8
        d.hkl.stored_string[] = "0 0 1"
        d.bfile.stored_string[] = joinpath(STRUCTS, "GaAs_conv.vasp")
        d.ab.stored_string[] = "1 2"
        d.build_sl.clicks[] += 1
        @test formula(v.crystal) == "Si8Ga8As8"
        d.strain.active[] = false
        undo!(v)
        d.build_sl.clicks[] += 1                                         # 4 % mismatch, no strain
        @test natoms(v.crystal) == 8 && startswith(v.status[], "Could not build")
    end

    @testset "docked tools" begin
        v = visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false)
        nblocks, nboxes = length(v.figure.content), length(v.textboxes)
        @test haskey(v.dialogs, :panel) && open_tool(v) === nothing
        d = AtomBond.supercell_dialog(v)                                     # opens on the right, in the main figure
        @test open_tool(v) === :supercell && !haskey(v.dialogs, :panel)
        @test any(t -> t === d.repeats, v.textboxes)                # shortcuts pause while typing in it
        @test AtomBond.supercell_dialog(v) === nothing && occursin("already open", v.status[])
        AtomBond.shift_dialog(v)                                             # another tool replaces it
        @test open_tool(v) === :shift && !any(t -> t === d.repeats, v.textboxes)
        @test close_tool!(v) && open_tool(v) === nothing && haskey(v.dialogs, :panel)
        @test length(v.figure.content) == nblocks && length(v.textboxes) == nboxes   # nothing left behind
        @test !close_tool!(v)
        # Long lines of a tool wrap at the column width.
        @test AtomBond.wrap_text("aa bb cc dd", 5) == "aa bb\ncc dd"
        @test AtomBond.wrap_text("x\nlongword y", 4) == "x\nlongword\ny"
        AtomBond.library_dialog(v)
        info = only(filter(b -> b isa AtomBond.Makie.Label && occursin("DOI", b.text[]), v.figure.content))
        @test occursin('\n', info.text[]) && all(l -> length(l) < 80, split(info.text[], '\n'))
        close_tool!(v)
        # WSL: software rendering unless asked for the GPU or a driver is chosen.
        @test !AtomBond.wsl_software_rendering!(Dict("ATOMBOND_GPU" => "1"))
        @test !AtomBond.wsl_software_rendering!(Dict("LIBGL_ALWAYS_SOFTWARE" => "0"))
        env = Dict{String,String}()
        wsl = Sys.islinux() && occursin("microsoft", lowercase(read("/proc/version", String)))
        @test AtomBond.wsl_software_rendering!(env) == wsl && haskey(env, "LIBGL_ALWAYS_SOFTWARE") == wsl
    end

    @testset "tabs" begin
        si = joinpath(STRUCTS, "Si_conv.vasp")
        v = visualize(si; show=false)
        @test ntabs(v) == 1 && [b.label[] for b in v.dialogs[:tabbar].buttons] == ["Si_conv.vasp"]
        select!(v, 1); AtomBond.delete_selected!(v)                 # an edit in tab 1
        @test formula(v.crystal) == "Si7" && v.dialogs[:tabbar].buttons[1].label[] == "Si_conv.vasp *"
        @test open_library!(v, "Water") && ntabs(v) == 2 && active_tab(v) == 2 && formula(v.crystal) == "O1H2"
        @test isempty(v.history)                                    # each tab has its own undo history
        select!(v, 1); select!(v, 2; add=true)
        switch_tab!(v, 1)
        @test formula(v.crystal) == "Si7" && length(v.history) == 1 && isempty(v.selection)
        @test undo!(v) && formula(v.crystal) == "Si8"
        switch_tab!(v, 2)
        @test formula(v.crystal) == "O1H2" && length(v.selection) == 2   # its selection comes back
        # A file already open (and unchanged) is shown in its tab, not opened twice.
        @test open_structure!(v, si) && ntabs(v) == 2 && active_tab(v) == 1
        # A structure tool follows the tab and names it.
        AtomBond.supercell_dialog(v)
        @test open_tool(v) === :supercell
        switch_tab!(v, 2)
        @test open_tool(v) === :supercell && any(b -> b isa AtomBond.Makie.Label && b.text[] == "on Water",
                                                 v.figure.content)
        close_tool!(v)
        # Closing: unsaved changes need a second close; the last tab leaves an empty cell.
        select!(v, 3); AtomBond.delete_selected!(v)
        @test !close_tab!(v) && startswith(v.status[], "Unsaved changes in Water *")
        @test close_tab!(v) && ntabs(v) == 1 && formula(v.crystal) == "Si8"
        @test close_tab!(v) && ntabs(v) == 1 && natoms(v.crystal) == 0 && AtomBond.is_blank(v)
        @test open_library!(v, "GaN") && ntabs(v) == 1                   # replaces the blank tab
        # new_tab!() starts an empty cell with the New structure tool; Create fills that tab.
        new_tab!(v)
        @test ntabs(v) == 2 && open_tool(v) === :new && AtomBond.is_blank(v)
        @test AtomBond.create_structure!(v, "Diamond", ["C"], (3.567, 3.567, 3.567)) && ntabs(v) == 2 && formula(v.crystal) == "C8"
        # Exit warns about unsaved changes in any tab.
        switch_tab!(v, 1); edit!(c -> supercell(c, 1, 1, 2), v, "double c")
        switch_tab!(v, 2)
        @test !exit!(v) && occursin("GaN", v.status[])
        # Shortcuts: Ctrl+Tab next tab, Ctrl+W close.
        M = AtomBond.Makie
        ev = M.events(v.lscene.scene)
        key(k) = (ev.keyboardbutton[] = M.KeyEvent(k, M.Keyboard.press); ev.keyboardbutton[] = M.KeyEvent(k, M.Keyboard.release))
        ev.keyboardbutton[] = M.KeyEvent(M.Keyboard.left_control, M.Keyboard.press)
        key(M.Keyboard.tab)
        @test active_tab(v) == 1
        key(M.Keyboard.w)
        @test ntabs(v) == 2 && occursin("Unsaved changes", v.status[])     # modified: warns first
        ev.keyboardbutton[] = M.KeyEvent(M.Keyboard.left_control, M.Keyboard.release)
    end

    @testset "oxygen bridge" begin
        V = AtomBond.Vec3
        a, b, c = V(0, 0, 0), V(2, 0, 0), V(1, 3, 0)
        @test bridge_position(a, b, c; offset=0.5) ≈ V(1, 0.5, 0)              # towards c
        @test bridge_position(a, b, c; offset=0.5, flip=true) ≈ V(1, -0.5, 0)  # -v
        @test bridge_position(a, b, V(1, 3, 4); offset=1) ≈ V(1, 0.6, 0.8)        # 1 Å towards the third atom
        p = bridge_position(a, b, V(5, 2, 1); offset=0.3)
        @test abs(dot(p - (a + b) / 2, b - a)) < 1e-12                          # perpendicular to the bond
        @test abs(dot(p - a, cross(b - a, V(5, 2, 1) - a))) < 1e-12              # in the plane of the three atoms
        @test AtomBond.bridge_max_offset(a, b, c) ≈ 3 && AtomBond.bridge_max_offset(a, b, V(3, 1, 0)) ≈ 1 / 3
        @test_throws ArgumentError bridge_position(a, b, V(5, 0, 0))              # third atom on the bond line
        @test_throws ArgumentError bridge_position(a, b, c; offset=-1)
        # bridges in molecules: three-membered rings keep their bonds, a bridging H does not make a B–B bond
        @test length(find_bonds(load_library("Cyclopropane"))) == 9 && length(find_bonds(load_library("Ethylene oxide (oxirane)"))) == 7
        b2h6 = load_library("Diborane")
        @test count(x -> b2h6.species[x.i] == b2h6.species[x.j] == :B, find_bonds(b2h6)) == 0 &&
              count(x -> b2h6.species[x.i] == b2h6.species[x.j] == :B, find_bonds(b2h6; bridge_angle=nothing)) == 1

        # In the window: Si, bond centre–neighbour 1, side from neighbour 2
        v = visualize(load_library("Silicon"); show=false)
        S = AtomBond.Shift
        bs = [x for x in find_bonds(v.crystal) if x.shift == S(0, 0, 0)]
        centre = argmax(i -> count(x -> i in (x.i, x.j), bs), 1:natoms(v.crystal))
        nb = [x.i == centre ? x.j : x.i for x in bs if centre in (x.i, x.j)]
        @test "Oxygen bridge…" in [it.label for it in AtomBond.menubar_items(v, "Build")]
        d = AtomBond.bridge_dialog(v)
        @test startswith(d.preview[], "Select three atoms")
        d.add.clicks[] = 1
        @test startswith(v.status[], "Select three atoms") && natoms(v.crystal) == 8
        v.selection = [(centre, S(0, 0, 0)), (nb[1], S(0, 0, 0)), (nb[2], S(0, 0, 0))]; AtomBond.update_selection!(v)
        @test occursin("Si–O 1.278 Å", replace(d.preview[], r"Si\d+" => "Si")) && occursin("Inside the triangle", d.preview[])
        @test "Add oxygen bridge" in [it.label for it in AtomBond.menu_items(v)]
        r = cartesian_positions(v.crystal)
        d.add.clicks[] = 1
        @test formula(v.crystal) == "Si8O1" && startswith(v.status[], "Add O bridge in")
        # the bridged Si–Si bond is gone; Si–O–Si instead, every Si still has four bonds
        nbonds(c, i) = count(x -> i in (x.i, x.j), find_bonds(c))
        pair(x) = Set((x.i, x.j))
        si_o = natoms(v.crystal)
        @test !any(x -> pair(x) == Set((centre, nb[1])) && x.shift == S(0, 0, 0), find_bonds(v.crystal))
        @test nbonds(v.crystal, si_o) == 2 && nbonds(v.crystal, centre) == 4 && nbonds(v.crystal, nb[1]) == 4
        @test any(x -> pair(x) == Set((centre, nb[1])), find_bonds(v.crystal; bridge_angle=nothing))   # old rule
        @test length(v.display.bonds) == length(find_bonds(v.crystal)) == 16 + 2 - 1      # drawn the same way
        o = cartesian_positions(v.crystal)[end]
        mid = (r[centre] + r[nb[1]]) / 2
        @test norm(o - mid) ≈ 0.5 && dot(o - mid, r[nb[2]] - mid) > 0                # towards the third atom
        @test occursin("is already there", d.preview[])
        # flip: the other side, then back; each step undoable
        @test "Flip oxygen bridge" in [it.label for it in AtomBond.menu_items(v)]
        d.flip_last.clicks[] = 1
        o2 = cartesian_positions(v.crystal)[end]
        @test o2 ≈ 2mid - o && startswith(v.status[], "Flip the O bridge")
        undo!(v)
        @test cartesian_positions(v.crystal)[end] ≈ o
        undo!(v)
        @test formula(v.crystal) == "Si8"
        d.flip_last.clicks[] = 1
        @test startswith(v.status[], "Nothing to flip")
        # options: other side up front, another element, a larger offset outside the triangle
        three = [(centre, S(0, 0, 0)), (nb[1], S(0, 0, 0)), (nb[2], S(0, 0, 0))]
        v.selection = copy(three); AtomBond.update_selection!(v)          # undo clears the selection
        d.flip.active[] = true
        @test occursin("other side of the bond", d.preview[])
        d.flip.active[] = false
        AtomBond.settext!(d.offset, "2")
        @test occursin("Outside the triangle", d.preview[])
        AtomBond.settext!(d.offset, "0.6"); AtomBond.settext!(d.element, "N")
        d.add.clicks[] = 1
        @test formula(v.crystal) == "Si8N1" && norm(cartesian_positions(v.crystal)[end] - mid) ≈ 0.6
        undo!(v)
        v.selection = copy(three); AtomBond.update_selection!(v)
        # the right-click item uses the window's settings (and defaults without it)
        @test AtomBond.bridge_settings(v) == (; element=:N, offset=0.6, flip=false)
        @test AtomBond.add_bridge!(v; AtomBond.bridge_settings(v)...) && formula(v.crystal) == "Si8N1"
        # collinear third atom
        w = visualize(load_library("CO2"); show=false)
        w.selection = [(1, S(0, 0, 0)), (2, S(0, 0, 0)), (3, S(0, 0, 0))]; AtomBond.update_selection!(w)
        @test !AtomBond.add_bridge!(w) && occursin("on the line of the bond", w.status[])
    end

    @testset "export" begin
        si = load_structure(joinpath(STRUCTS, "Si_conv.vasp"))
        gan = load_structure(joinpath(STRUCTS, "GaN_wurtzite.vasp"))
        rut = load_structure(joinpath(STRUCTS, "TiO2_rutile.cif"))
        @test kpoint_grid(si) == (5, 5, 5)                       # |b| = 2π/5.431 = 1.157 Å⁻¹
        @test kpoint_grid(primitive_cell(si)) == (9, 9, 9)
        @test kpoint_grid(si; kspacing=0.5) == (3, 3, 3)
        @test kpoint_grid(slab(si, (1, 1, 1); layers=2, vacuum=15))[3] == 1   # vacuum
        @test kpoint_grid(supercell(si, 1, 1, 4)) == (5, 5, 2)                # bulk, no vacuum
        @test kpoint_grid(load_structure(joinpath(STRUCTS, "water.xyz"))) == (1, 1, 1)
        @test atomic_mass(:Si) ≈ 28.085
        @test_throws ArgumentError kpoint_grid(si; kspacing=0)

        # CIF (P1) round trip, including a slab.
        for c in (gan, rut, slab(si, (1, 1, 0); layers=2, vacuum=10))
            io = IOBuffer(); write_cif(io, c)
            back = read_cif(IOBuffer(String(take!(io))))
            @test back.species == c.species
            @test all(isapprox.(back.frac, c.frac; atol=1e-7))
            @test all(isapprox.(lattice_parameters(back.lattice), lattice_parameters(c.lattice); atol=1e-5))
        end

        # Quantum ESPRESSO: parse the namelist values and cards back.
        io = IOBuffer(); write_qe(io, gan; ecutwfc=50, pseudo="{El}.pbe.UPF")
        qe = String(take!(io))
        @test occursin("nat = 4", qe) && occursin("ntyp = 2", qe) && occursin("ecutwfc = 50.0", qe)
        @test occursin(r"Ga\s+69\.7230\s+Ga\.pbe\.UPF", qe) && occursin(r"N\s+14\.0070\s+N\.pbe\.UPF", qe)
        lines = split(qe, '\n')
        k = findfirst(==("CELL_PARAMETERS angstrom"), lines)
        cell = [parse.(Float64, split(lines[k+i])) for i in 1:3]
        @test hcat(cell...) ≈ gan.lattice.matrix
        k = findfirst(==("ATOMIC_POSITIONS crystal"), lines)
        @test [split(lines[k+i])[1] for i in 1:4] == ["Ga", "Ga", "N", "N"]
        @test parse.(Float64, split(lines[k+3])[2:4]) ≈ gan.frac[3]
        k = findfirst(==("K_POINTS automatic"), lines)
        @test split(lines[k+1]) == string.([kpoint_grid(gan)..., 0, 0, 0])

        # CASTEP .cell
        io = IOBuffer(); write_castep(io, rut)
        ce = String(take!(io))
        lines = split(ce, '\n')
        k = findfirst(==("%BLOCK LATTICE_CART"), lines)
        @test hcat([parse.(Float64, split(lines[k+1+i])) for i in 1:3]...) ≈ rut.lattice.matrix
        k = findfirst(==("%BLOCK POSITIONS_FRAC"), lines)
        @test count(l -> startswith(strip(l), "Ti "), lines[k+1:k+6]) == 2
        @test occursin("KPOINTS_MP_GRID " * join(kpoint_grid(rut), " "), ce)

        # export_inputs for every format, and save_structure by extension.
        mktempdir() do dir
            p = export_inputs(gan, :vasp, dir, "GaN")
            @test basename.(p) == ["POSCAR", "KPOINTS", "POTCAR.order"]
            @test read_poscar(p[1]).species == gan.species
            @test strip(readlines(p[2])[4]) == join(kpoint_grid(gan), " ")
            @test strip(read(p[3], String)) == "Ga N"
            for (fmt, ext) in [(:qe, ".in"), (:abinit, ".abi"), (:siesta, ".fdf"), (:castep, ".cell"), (:cif, ".cif"),
                               (:xyz, ".xyz")]
                @test isfile(only(export_inputs(gan, fmt, dir, "GaN")))
                @test endswith(only(export_inputs(gan, fmt, dir, "GaN")), ext)
            end
            @test occursin("pseudos \"Ga.psp8, N.psp8\"", read(joinpath(dir, "GaN.abi"), String))   # ABINIT default
            @test_throws ArgumentError export_inputs(gan, :wien2k, dir, "GaN")
            # ABINIT and SIESTA inputs read back as the same structure
            for (ext, writer) in ((".abi", write_abinit), (".fdf", write_siesta))
                for x in (gan, rut, si)
                    writer(joinpath(dir, "rt" * ext), x)
                    back = load_structure(joinpath(dir, "rt" * ext))
                    @test (ext, back.species) == (ext, x.species)
                    @test back.lattice.matrix ≈ x.lattice.matrix atol=1e-8
                    @test all(isapprox.(back.frac, x.frac; atol=1e-8))
                end
            end
            ab = sprint(io -> write_abinit(io, gan; kspacing=0.5, ecutwfc=50))
            @test occursin("ngkpt " * join(kpoint_grid(gan; kspacing=0.5), " "), ab) && occursin("ecut 25.00", ab) &&
                  occursin("znucl 31 7", ab) && occursin("typat 1 1 2 2", ab) && !occursin("ionmov", ab)
            fd = sprint(io -> write_siesta(io, gan))
            @test occursin("  1 31 Ga\n  2 7 N", fd) && occursin("%block kgrid_Monkhorst_Pack", fd) && !occursin("MD.TypeOfRun", fd)
            # constraints and spin
            gc = with_constraints(fix_atoms(gan, [1, 2]); cell=(1, 1, 1, 0, 0, 0))
            ab = sprint(io -> write_abinit(io, gc))
            @test occursin("ionmov 2", ab) && occursin("optcell 1", ab) && occursin("natfix 2\niatfix 1 2", ab)
            fd = sprint(io -> write_siesta(io, gc))
            @test occursin("MD.TypeOfRun CG", fd) && occursin("MD.VariableCell T", fd) &&
                  occursin("atom 1 2", fd) && occursin("cell constraints are not applied", fd)
            mag = with_magmoms(si, fill(1.0, natoms(si)))
            @test occursin("nsppol 2\nspinat", sprint(io -> write_abinit(io, mag)))
            @test occursin("Spin polarized\n%block DM.InitSpin\n  1 1.0000", sprint(io -> write_siesta(io, mag)))
            # DFTK.jl script: valid Julia, the same atoms, the k-grid, cutoff and pseudopotentials
            dj = sprint(io -> write_dftk(io, gan; kspacing=0.5, ecutwfc=50))
            @test Meta.parseall(dj).head === :toplevel && !any(e -> e isa Expr && e.head in (:error, :incomplete),
                                                              Meta.parseall(dj).args)
            rows = [m.captures for m in eachmatch(r"Atom\(:(\w+), \[([-\d.]+), ([-\d.]+), ([-\d.]+)\]", dj)]
            @test Symbol.(first.(rows)) == gan.species
            @test all(isapprox.([parse.(Float64, r[2:4]) for r in rows], cartesian_positions(gan); atol=1e-8))
            @test occursin("kgrid=[" * join(kpoint_grid(gan; kspacing=0.5), ", ") * "]", dj) && occursin("Ecut=25.00", dj)
            @test occursin("PseudoFamily(\"$(AtomBond.DFTK_PSEUDO_FAMILY)\")", dj) && occursin("self_consistent_field", dj)
            @test !occursin("minimize_energy!", dj) && !occursin("magnetic_moment", dj)
            @test occursin("Dict(:Ga => \"Ga.upf\", :N => \"N.upf\")", sprint(io -> write_dftk(io, gan; pseudo="{El}.upf")))
            dj = sprint(io -> write_dftk(io, gc))
            @test occursin("using GeometryOptimization", dj) && occursin("variablecell=true", dj) &&
                  count("clamp=true", dj) == 2 && occursin("cell constraints are not applied", dj) &&
                  occursin("symmetries=false", dj)
            @test occursin("variablecell=false", sprint(io -> write_dftk(io, with_constraints(gc; cell=(0, 0, 0, 0, 0, 0)))))
            @test count("magnetic_moment=1.0000", sprint(io -> write_dftk(io, mag))) == natoms(si)
            p = only(export_inputs(gan, :dftk, dir, "GaN"))                   # default .UPF pattern → PseudoDojo family
            @test endswith(p, "GaN.jl") && occursin("PseudoFamily(", read(p, String))
            @test_throws ArgumentError export_inputs(gan, :qe, dir, "GaN"; pseudo=AtomBond.DFTK_PSEUDO_FAMILY)
            @test natoms(load_structure(save_structure(joinpath(dir, "s.fdf"), si))) == natoms(si)
            @test natoms(load_structure(save_structure(joinpath(dir, "s.abi"), si))) == natoms(si)
            @test natoms(load_structure(save_structure(joinpath(dir, "rt.cif"), rut))) == 6
            @test isfile(save_structure(joinpath(dir, "x.in"), si; kspacing=0.4))
        end
    end

    @testset "AtomsBase" begin
        gan = load_structure(joinpath(STRUCTS, "GaN_wurtzite.vasp"))
        sys = to_atomsbase(gan)
        @test length(sys) == 4
        @test AtomBond.AtomsBase.atomic_symbol(sys, 3) == :N
        back = from_atomsbase(sys)
        @test back.species == gan.species && back.pbc == gan.pbc
        @test all(isapprox.(back.frac, gan.frac; atol=1e-10))
        @test back.lattice.matrix ≈ gan.lattice.matrix
        w = to_atomsbase(load_structure(joinpath(STRUCTS, "water.xyz")))
        @test all(.!AtomBond.AtomsBase.periodicity(w))
        # magnetic moments and fully fixed atoms travel as the atom properties DFTK and GeometryOptimization use
        mc = with_magmoms(fix_atoms(gan, [2]), [1.5, NaN, 0.0, NaN])
        ms = to_atomsbase(mc)
        @test ms[1][:magnetic_moment] == 1.5 && !haskey(ms[2], :magnetic_moment) && ms[2][:clamp] === true
        back = from_atomsbase(ms)
        @test isequal(back.magmoms, mc.magmoms) && back.fixed == mc.fixed
        @test !has_magmoms(from_atomsbase(sys)) && !has_constraints(from_atomsbase(sys))
    end

    @testset "export dialog" begin
        v = visualize(joinpath(STRUCTS, "GaN_wurtzite.vasp"); show=false, panel=false)
        d = AtomBond.export_dialog(v)
        @test d.mesh[] == "Mesh for this structure: 10 × 10 × 5"
        d.kspacing.stored_string[] = "0.5"
        @test d.mesh[] == "Mesh for this structure: 5 × 5 × 3"
        mktempdir() do dir
            d.folder.stored_string[] = dir
            d.name.stored_string[] = "gan"
            @test d.format.selection[] === :qe                 # open-source codes first
            d.format.i_selected[] = 2                         # ABINIT: the pseudopotential pattern follows
            @test AtomBond.shown_text(d.pseudo) == "{El}.psp8"
            d.go.clicks[] += 1
            @test isfile(joinpath(dir, "gan.abi")) && occursin("ngkpt 5 5 3", read(joinpath(dir, "gan.abi"), String))
            d.format.i_selected[] = findfirst(o -> o[2] === :dftk, AtomBond.EXPORT_FORMATS)   # DFTK: a pseudo family
            @test AtomBond.shown_text(d.pseudo) == AtomBond.DFTK_PSEUDO_FAMILY
            d.go.clicks[] += 1
            @test isfile(joinpath(dir, "gan.jl")) && occursin("kgrid=[5, 5, 3]", read(joinpath(dir, "gan.jl"), String))
            d.format.i_selected[] = 1                         # Quantum ESPRESSO
            @test AtomBond.shown_text(d.pseudo) == "{El}.UPF"
            d.go.clicks[] += 1
            @test isfile(joinpath(dir, "gan.in"))
            @test occursin("K_POINTS automatic\n  5 5 3 0 0 0", read(joinpath(dir, "gan.in"), String))
            @test startswith(v.status[], "Wrote gan.in")
        end
        # The mesh follows structure changes in the main window.
        AtomBond.edit!(c -> supercell(c, 2, 2, 1), v, "make supercell")
        @test d.mesh[] == "Mesh for this structure: 3 × 3 × 3"
    end

    @testset "load_structure" begin
        # (other files, e.g. images saved from the app, are skipped)
        for f in filter(f -> isfile(f) && structure_format(f) !== nothing, readdir(STRUCTS; join=true))
            @test natoms(load_structure(f)) > 0
        end
        @test_throws ArgumentError load_structure(joinpath(STRUCTS, "nope.cif"))
        @test structure_format("CONTCAR_relaxed") == :poscar && structure_format("x.CELL") == :castep_cell
        @test structure_format("run.geom") == :castep_geom && structure_format("notes.txt") === nothing
    end

    @testset "CASTEP readers" begin
        # .cell: LATTICE_ABC, lower-case blocks, comments, labelled species
        c = load_structure(joinpath(STRUCTS, "GaAs_abc.cell"))
        @test c.species == [:Ga, :As] && c.frac[2] ≈ [0.25, 0.25, 0.25]
        @test spacegroup(c) == (216, "F-43m")
        @test all(isapprox.(lattice_parameters(c.lattice), (3.997, 3.997, 3.997, 60, 60, 60); atol=1e-8))
        # round trip through our writer (LATTICE_CART, POSITIONS_FRAC)
        for x in (load_structure(joinpath(STRUCTS, "GaN_wurtzite.vasp")), load_structure(joinpath(STRUCTS, "TiO2_rutile.cif")))
            io = IOBuffer(); write_castep(io, x)
            back = read_castep_cell(IOBuffer(String(take!(io))))
            @test back.species == x.species[AtomBond.grouped_order(x)]
            @test all(isapprox.(back.frac, x.frac[AtomBond.grouped_order(x)]; atol=1e-9))
            @test back.lattice.matrix ≈ x.lattice.matrix
        end
        # POSITIONS_ABS in Bohr, species given by atomic number
        b = read_castep_cell(IOBuffer("""
            %BLOCK LATTICE_CART
            bohr
              10.2631 0 0
              0 10.2631 0
              0 0 10.2631
            %ENDBLOCK LATTICE_CART
            %block positions_abs
            bohr
              14 0 0 0
              Si 2.565775 2.565775 2.565775
            %endblock positions_abs
            """))
        @test b.species == [:Si, :Si] && b.frac[2] ≈ [0.25, 0.25, 0.25]
        @test lattice_parameters(b.lattice)[1] ≈ 10.2631 * AtomBond.BOHR
        @test_throws ArgumentError read_castep_cell(IOBuffer("%block positions_frac\nSi 0 0 0\n%endblock positions_frac\n"))

        # .geom: frames in Bohr; load_structure takes the last (relaxed) one
        frames = read_castep_geom_frames(joinpath(STRUCTS, "Si_relax.geom"))
        @test length(frames) == 2
        @test frames[1].frac[2] ≈ fill(2.4 / (4 * 2.56557), 3) atol = 1e-6
        @test frames[2].frac[2] ≈ [0.25, 0.25, 0.25]
        @test lattice_parameters(frames[2].lattice)[1] ≈ 5.13114 * sqrt(2) * AtomBond.BOHR
        @test occursin("step 1", frames[2].title) && occursin("E = -391.34", frames[2].title)
        @test load_structure(joinpath(STRUCTS, "Si_relax.geom")).frac == frames[2].frac
        @test read_castep_geom(joinpath(STRUCTS, "Si_relax.geom"); frame=1).frac == frames[1].frac
        @test spacegroup(frames[2]) == (227, "Fd-3m")
    end

    @testset "edit atoms" begin
        gaas = load_structure(joinpath(STRUCTS, "GaAs_conv.vasp"))
        # set_position: fractional and Cartesian, overlap refused
        c = set_position(gaas, 5, (0.3, 0.25, 0.25))
        @test c.frac[5] ≈ [0.3, 0.25, 0.25] && c.frac[1:4] == gaas.frac[1:4]
        c = set_position(gaas, 5, (1.5, 1.5, 1.5); cartesian=true)
        @test cartesian_positions(c)[5] ≈ [1.5, 1.5, 1.5]
        @test_throws ArgumentError set_position(gaas, 5, (0, 0, 0))     # on top of Ga1
        @test_throws ArgumentError set_position(gaas, 9, (0, 0, 0))
        # replace_species
        p = replace_species(gaas, "As", "p")
        @test formula(p) == "Ga4P4" && p.frac == gaas.frac
        @test_throws ArgumentError replace_species(gaas, :Si, :Ge)       # no Si
        @test_throws ArgumentError replace_species(gaas, :Ga, :Ga)
        @test_throws ArgumentError replace_species(gaas, :Ga, :Xx)

        # Window
        v = visualize(joinpath(STRUCTS, "GaAs_conv.vasp"); show=false)
        d = AtomBond.atoms_dialog(v)
        @test startswith(d.who[], "None") && d.from.options[] == ["Ga", "As"] && d.count_text[] == "4 atoms of Ga"
        select!(v, 5)
        @test d.who[] == "As5" && d.element.displayed_string[] == "As" && d.from.selection[] == "As"
        @test d.coords.displayed_string[] == "0.250000  0.250000  0.250000"
        d.mode.i_selected[] = 2
        @test d.coords.displayed_string[] == "1.413250  1.413250  1.413250"
        # typed (Enter not pressed) values are used; one undoable edit
        d.element.displayed_string[] = "in"
        d.coords.displayed_string[] = "1.5 1.5 1.5"
        d.apply.clicks[] += 1
        @test v.crystal.species[5] == :In && cartesian_positions(v.crystal)[5] ≈ [1.5, 1.5, 1.5]
        @test v.selection == [(5, AtomBond.Shift(0, 0, 0))] && length(v.history) == 1
        d.apply.clicks[] += 1
        @test v.status[] == "Nothing to change"
        undo!(v)
        @test v.crystal.species[5] == :As && v.crystal.frac[5] ≈ [0.25, 0.25, 0.25]
        # an image copy moves to the typed position of that copy
        d.mode.i_selected[] = 1
        clear_selection!(v); select!(v, 1, AtomBond.Shift(1, 0, 0))
        @test endswith(d.who[], "(copy in a neighbouring cell)") && startswith(d.coords.displayed_string[], "1.000000")
        d.coords.displayed_string[] = "1 0 0.1"
        d.apply.clicks[] += 1
        @test v.crystal.frac[1] ≈ [0, 0, 0.1]
        d.coords.displayed_string[] = "0.25 0.25 0.25"                   # onto As5
        d.apply.clicks[] += 1
        @test startswith(v.status[], "Could not edit Ga1") && v.crystal.frac[1] ≈ [0, 0, 0.1]
        # several atoms: element for all
        select!(v, 2; add=true); select!(v, 6; add=true)
        @test d.who[] == "3 atoms: Ga×2, As×1" && d.element.displayed_string[] == ""
        d.apply.clicks[] += 1
        @test v.status[] == "Enter an element for the selected atoms"
        d.element.displayed_string[] = "Al"
        d.apply.clicks[] += 1
        @test v.crystal.species[[1, 2, 6]] == [:Al, :Al, :Al]
        # replace a species everywhere
        d.from.i_selected[] = findfirst(==("As"), d.from.options[])
        d.to.displayed_string[] = "Sb"
        d.replace_all.clicks[] += 1
        @test formula(v.crystal) == "Al3Ga2Sb3" && startswith(v.status[], "Replace all As with Sb (3 atoms)")
        @test d.from.options[] == ["Al", "Ga", "Sb"]
        d.to.displayed_string[] = "Qq"
        d.replace_all.clicks[] += 1
        @test v.status[] == "Unknown element 'Qq'"
    end

    @testset "bonds kept after substitution" begin
        gaas = load_structure(joinpath(STRUCTS, "GaAs_conv.vasp"))
        gac = replace_species(gaas, :As, :C)
        @test isempty(find_bonds(gac))                    # Ga–C at the GaAs distance (2.45 Å) > 2.28 Å
        cut = Dict{Tuple{Symbol,Symbol},Float64}()
        @test AtomBond.inherit_bond_cutoffs!(cut, gaas, gac) == [(:C, :Ga)]
        @test cut[(:C, :Ga)] ≈ bond_cutoff(:Ga, :As) ≈ 1.15 * (1.22 + 1.19)
        @test length(find_bonds(gac; cutoffs=cut)) == 16
        @test length(find_bonds(supercell(gac, 2, 2, 2); cutoffs=cut)) == 128
        @test isempty(AtomBond.inherit_bond_cutoffs!(cut, gaas, replace_species(gaas, :As, :Sb)))   # Ga–Sb is long anyway
        # The viewer does this on any element-only edit and says so.
        v = visualize(joinpath(STRUCTS, "GaAs_conv.vasp"); show=false)
        edit!(c -> replace_species(c, :As, :C), v, "replace As with C")
        @test length(v.display.bonds) == 16 && occursin("bonds kept: C–Ga ≤ 2.77 Å", v.status[])
        select!(v, 1)
        @test endswith(v.info[], "4 bonds")
        select!(v, 5)                                      # single-atom path (panel Replace)
        substitute_selected = AtomBond.substitute_selected!
        substitute_selected(v, :N)
        @test length(v.display.bonds) == 16 && haskey(v.settings.bond_cutoffs, (:Ga, :N))
        # moves don't add cutoffs
        n = length(v.settings.bond_cutoffs)
        edit!(c -> move_atoms(c, [5], (0.1, 0, 0)), v, "move")
        @test length(v.settings.bond_cutoffs) == n
        # other pairs unaffected: Si unchanged
        si = load_structure(joinpath(STRUCTS, "Si_conv.vasp"))
        @test length(find_bonds(si; cutoffs=v.settings.bond_cutoffs)) == 16
    end

    @testset "panel text boxes" begin
        M = AtomBond.Makie
        @test AtomBond.tryparse3(Int, "2x2x1") == [2, 2, 1] && AtomBond.tryparse3(Int, " 1, 1 ,3 ") == [1, 1, 3]
        @test AtomBond.tryparse3(Int, "2×3×4") == [2, 3, 4] && AtomBond.tryparse3(Int, "2 2") === nothing
        @test AtomBond.tryparse3(Float64, "0.1, -2e-1; 3") == [0.1, -0.2, 3.0]
        @test AtomBond.tryparse3(Float64, "1x2x3") === nothing
        v = visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false)
        box(p) = only(filter(t -> t.placeholder[] == p, v.textboxes))
        button(l) = only(filter(b -> b isa M.Button && b.label[] == l, v.figure.content))
        sc, go = box("n1 n2 n3"), button("Supercell")
        go.clicks[] += 1                                   # empty box: a red hint
        @test natoms(v.crystal) == 8 && v.status[] == "Type the repeats n1 n2 n3 first"
        @test AtomBond.status_color(v.status[]) == RGBf(0.75, 0.1, 0.1)
        sc.displayed_string[] = "2x2x1"                    # typed, Enter not pressed
        go.clicks[] += 1
        @test natoms(v.crystal) == 32 && AtomBond.status_color(v.status[]) == RGBf(0.3, 0.3, 0.3)
        sc.displayed_string[] = "2 2"
        go.clicks[] += 1
        @test natoms(v.crystal) == 32 && startswith(v.status[], "Cannot use '2 2'")
        select!(v, 1)
        box("move by dx dy dz (Å)").displayed_string[] = "0.1, 0, 0"
        button("Move").clicks[] += 1
        @test startswith(v.status[], "Move 1 atom")
        # dialogs read typed text too
        d = AtomBond.surface_dialog(v)
        undo!(v); undo!(v)
        d.hkl.displayed_string[] = "1,1,1"
        d.layers.displayed_string[] = "2"
        d.build_slab.clicks[] += 1
        @test startswith(v.status[], "Cleave (1 1 1) slab, 2 layers")
    end

    @testset "file browser" begin
        M = AtomBond.Makie
        @test [e[1] for e in AtomBond.browser_entries(joinpath(STRUCTS, ".."))] == ["structures"]
        tmp = mktempdir()
        mkdir(joinpath(tmp, "b_dir")); mkdir(joinpath(tmp, "A_dir")); mkdir(joinpath(tmp, ".hidden"))
        for f in ("notes.txt", "POSCAR", "z.cif", "a.geom")
            touch(joinpath(tmp, f))
        end
        @test AtomBond.browser_entries(tmp) ==
              [("A_dir", true), ("b_dir", true), ("a.geom", false), ("POSCAR", false), ("z.cif", false)]
        @test length(AtomBond.browser_entries(tmp; all_files=true)) == 6

        v = visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false)
        @test v.folder == rstrip(normpath(STRUCTS), '/')
        b = AtomBond.open_dialog(v)
        @test b.dir == rstrip(normpath(STRUCTS), '/') && length(b.entries) == count(f -> isdir(joinpath(STRUCTS, f)) || structure_format(f) !== nothing, readdir(STRUCTS))
        AtomBond.GLMakie.activate!(ssao=false)
        M.colorbuffer(b.v.figure)                                  # lay out the window
        ev = M.events(b.v.figure.scene)
        press() = (ev.mousebutton[] = M.MouseButtonEvent(M.Mouse.left, M.Mouse.press);
                   ev.mousebutton[] = M.MouseButtonEvent(M.Mouse.left, M.Mouse.release))
        i = findfirst(e -> e[1] == "TiO2_rutile.cif", b.entries)
        ev.mouseposition[] = AtomBond.row_pixel(b, i)
        press()
        @test b.selected == i && startswith(b.preview[], "Ti2O4, 6 atoms")
        @test formula(v.crystal) == "Si8"                         # one click only selects
        b.last_click = (0, 0.0); press(); press()                 # double-click opens
        @test formula(v.crystal) == "Ti2O4" && v.original === v.crystal
        @test startswith(v.status[], "Opened TiO2_rutile.cif") && ntabs(v) == 2 && active_tab(v) == 2   # a new tab
        switch_tab!(v, 1)
        @test formula(v.crystal) == "Si8"
        # Keys and folders
        ev.keyboardbutton[] = M.KeyEvent(M.Keyboard.down, M.Keyboard.press)
        @test b.selected == i + 1
        ev.keyboardbutton[] = M.KeyEvent(M.Keyboard.backspace, M.Keyboard.press)
        @test b.dir == dirname(normpath(STRUCTS))
        @test AtomBond.navigate!(b, tmp) && length(b.entries) == 5
        AtomBond.select_entry!(b, 1)
        @test AtomBond.activate!(b) && b.dir == joinpath(tmp, "A_dir")    # enter a folder
        @test b.position[] == "No folders or structure files here"
        @test !AtomBond.navigate!(b, joinpath(tmp, "missing")) && startswith(v.status[], "Cannot open folder")
        AtomBond.navigate!(b, tmp)
        AtomBond.select_entry!(b, 4)                              # empty POSCAR
        @test startswith(b.preview[], "Cannot read")
        @test !AtomBond.activate!(b) && formula(v.crystal) == "Si8"
        # Scrolling a long list
        for k in 1:40; touch(joinpath(tmp, "f$k.xyz")); end
        AtomBond.navigate!(b, tmp)
        ev.mouseposition[] = AtomBond.row_pixel(b, 2)
        ev.scroll[] = (0.0, -1.0)
        @test b.offset == 3 && startswith(b.position[], "4–19 of 45")
        AtomBond.select_entry!(b, 45)
        @test b.offset == 45 - AtomBond.BROWSER_ROWS
        # folder text box
        b.folder.stored_string[] = STRUCTS
        @test b.dir == rstrip(normpath(STRUCTS), '/')
    end

    @testset "right-click menu" begin
        M = AtomBond.Makie
        v = visualize(joinpath(STRUCTS, "GaAs_conv.vasp"); show=false)
        AtomBond.GLMakie.activate!(ssao=false)
        M.colorbuffer(v.figure)
        ev = M.events(v.figure.scene)
        click(b) = (ev.mousebutton[] = M.MouseButtonEvent(b, M.Mouse.press);
                    ev.mousebutton[] = M.MouseButtonEvent(b, M.Mouse.release))
        labels() = [it.label for it in v.dialogs[:menu].items]
        # Nothing selected, nothing under the mouse: no menu (right-drag pans).
        ev.mouseposition[] = (5.0, 5.0)
        click(M.Mouse.right)
        @test !AtomBond.menu_isopen(v)
        # Right-click on an atom selects it and opens the menu.
        sc = v.lscene.scene
        atoms = AtomBond.drawn_atoms(v.display, v.settings)
        px(a) = M.project(sc, M.Point3f(a.position)) .+ M.origin(sc.viewport[])
        k = findfirst(a -> iszero(a.image) && (h = AtomBond.atom_at(v, px(a))) !== nothing && h.index == a.index, atoms)
        a = atoms[k]
        ev.mouseposition[] = Tuple(px(a))
        click(M.Mouse.right)
        @test v.selection == [(a.index, a.image)] && AtomBond.menu_isopen(v)
        el = v.crystal.species[a.index]
        other = only(filter(!=(el), species_list(v.crystal)))
        @test labels()[2:end] == ["Edit atoms…", "Shift atoms…", "Rotate atoms…", "Change to $other", "Select all $el", "Fix 1 atom", "Delete 1 atom",
                                  "Clear selection"]
        @test startswith(labels()[1], AtomBond.atom_label(v.crystal, a.index))
        # A click outside closes it and keeps the selection.
        ev.mouseposition[] = (5.0, 5.0)
        click(M.Mouse.left)
        @test !AtomBond.menu_isopen(v) && length(v.selection) == 1
        # Esc closes it and keeps the selection.
        AtomBond.open_context_menu!(v, (300, 500))
        ev.keyboardbutton[] = M.KeyEvent(M.Keyboard.escape, M.Keyboard.press)
        @test !AtomBond.menu_isopen(v) && length(v.selection) == 1
        # Clicking the title does nothing; clicking an item runs it.
        AtomBond.open_context_menu!(v, (300, 500))
        ev.mouseposition[] = AtomBond.menu_row_pixel(v, 1)
        click(M.Mouse.left)
        @test AtomBond.menu_isopen(v)
        ev.mouseposition[] = AtomBond.menu_row_pixel(v, 6)          # Select all
        @test v.dialogs[:menu].hover[] == 6
        click(M.Mouse.left)
        @test !AtomBond.menu_isopen(v) && length(AtomBond.selected_indices(v)) == count(==(el), v.crystal.species)
        AtomBond.open_context_menu!(v, (300, 500))
        @test labels()[1] == "4 atoms: $(el)×4" && "Delete 4 atoms" in labels()
        ev.mouseposition[] = AtomBond.menu_row_pixel(v, 5)          # Change to the other species
        click(M.Mouse.left)
        @test formula(v.crystal) == "$(other)8" && occursin("bonds kept", v.status[])
        undo!(v)
        # A right-click on an atom that is already selected keeps the selection.
        AtomBond.select_species!(v, el)
        ev.mouseposition[] = Tuple(px(a))
        click(M.Mouse.right)
        @test length(v.selection) == 4 && AtomBond.menu_isopen(v)
        clear_selection!(v)                                         # selection change closes it
        @test !AtomBond.menu_isopen(v)
        @test AtomBond.open_context_menu!(v, (300, 500)) === nothing         # nothing selected
    end

    @testset "prototypes" begin
        c = prototype_crystal("Diamond", ["Si"], (5.431, 5.431, 5.431))
        @test formula(c) == "Si8" && spacegroup(c)[1] == 227
        @test spacegroup(prototype_crystal("Zincblende", ["Ga", "As"], (5.65, 5.65, 5.65)))[1] == 216
        @test spacegroup(prototype_crystal("Rocksalt", ["Na", "Cl"], (5.64, 5.64, 5.64)))[1] == 225
        @test spacegroup(prototype_crystal("Fluorite", ["Ca", "F"], (5.46, 5.46, 5.46)))[1] == 225
        @test spacegroup(prototype_crystal("Perovskite ABX₃", ["Sr", "Ti", "O"], (3.9, 3.9, 3.9)))[1] == 221
        @test spacegroup(prototype_crystal("Hexagonal close-packed", ["Mg"], (3.21, 3.21, 5.21)))[1] == 194
        w = prototype_crystal("Wurtzite", ["Ga", "N"], (3.189, 3.189, 5.185))
        @test spacegroup(w)[1] == 186 && 1.9 < minimum(b -> b.length, find_bonds(w)) < 2.0
        @test all(==(4), coordination(w, find_bonds(w)))
        e = prototype_crystal("Empty cell", String[], (4, 5, 6), (90, 100, 90))
        @test natoms(e) == 0 && lattice_parameters(e.lattice)[5] ≈ 100
        @test_throws ArgumentError prototype_crystal("Rocksalt", ["Na"], (5.6, 5.6, 5.6))
        @test_throws ArgumentError prototype_crystal("Quasicrystal", ["Al"], (1, 1, 1))
        @test AtomBond.parse_abc("5.4") == (5.4, 5.4, 5.4) && AtomBond.parse_abc("3 5") == (3.0, 3.0, 5.0)
        @test AtomBond.parse_abc("1 2 3 4") === nothing && AtomBond.parse_abc("") === nothing
    end

    @testset "menu bar, new, save" begin
        M = AtomBond.Makie
        src = joinpath(STRUCTS, "GaAs_conv.vasp")
        dir = mktempdir()
        cp(src, joinpath(dir, "GaAs.vasp")); cp(joinpath(STRUCTS, "Si_Fd-3m.cif"), joinpath(dir, "Si.cif"))
        v = visualize(joinpath(dir, "GaAs.vasp"); show=false)
        @test v.document[] == "GaAs.vasp"
        # Menu bar: click a title, hover switches menus, click its title again closes.
        AtomBond.GLMakie.activate!(ssao=false)
        M.colorbuffer(v.figure)
        ev = M.events(v.figure.scene)
        click() = (ev.mousebutton[] = M.MouseButtonEvent(M.Mouse.left, M.Mouse.press);
                   ev.mousebutton[] = M.MouseButtonEvent(M.Mouse.left, M.Mouse.release))
        center(b) = (bb = b.layoutobservables.computedbbox[]; Tuple(M.origin(bb) .+ M.widths(bb) ./ 2))
        bar = v.dialogs[:menubar]
        ev.mouseposition[] = center(bar[1]); click()
        m = v.dialogs[:menu]
        @test AtomBond.menu_isopen(v) && m.owner == 1 && [it.label for it in m.items][1:2] == ["New structure…", "Open…"]
        ev.mouseposition[] = center(bar[2])
        @test m.owner == 2 && m.items[1].kind === :disabled           # Edit: nothing to undo
        click()
        @test !AtomBond.menu_isopen(v)
        # Hover highlights every item, also after shorter menus were open
        # (listeners of earlier menus used to break rows past their end).
        for t in (4, 3, 2, 1)
            AtomBond.open_menubar_menu!(v, t)
            for (j, it) in enumerate(m.items)
                it.kind === :item || continue
                ev.mouseposition[] = AtomBond.menu_row_pixel(v, j)
                @test m.hover[] == j
            end
        end
        AtomBond.close_context_menu!(v)
        ev.mouseposition[] = center(bar[2]); click()
        k = findfirst(it -> it.label == "Select all", m.items)
        @test AtomBond.menu_row_at(v, AtomBond.menu_row_pixel(v, 1)) == 0  # disabled Undo
        ev.mouseposition[] = AtomBond.menu_row_pixel(v, k); click()
        @test length(v.selection) == 8 && !AtomBond.menu_isopen(v)
        # Edits mark the file modified; Save writes it back (POSCAR).
        AtomBond.delete_selected!(v)
        @test v.document[] == "GaAs.vasp — modified"
        @test AtomBond.menubar_items(v, "Edit")[1].label == "Undo delete 8 atoms"
        undo!(v)
        @test v.document[] == "GaAs.vasp"
        edit!(c -> AtomBond.substitute(c, [1], :In), v, "change Ga1 to In")
        @test save!(v) === true && v.document[] == "GaAs.vasp"
        @test formula(load_structure(joinpath(dir, "GaAs.vasp"))) == "In1Ga3As4"
        @test AtomBond.recent_files()[1] == joinpath(dir, "GaAs.vasp")
        # A CIF from disk is not overwritten by Save: Save as opens instead.
        @test open_structure!(v, joinpath(dir, "Si.cif")) && v.document[] == "Si.cif"
        @test AtomBond.recent_files()[1] == joinpath(dir, "Si.cif")
        b = save!(v)
        @test b isa AtomBond.FileBrowser && b.dir == dir && AtomBond.shown_text(b.name) == "Si.cif"
        # Replacing needs a second click; the format menu sets the extension.
        @test !AtomBond.save_from_browser!(b) && startswith(v.status[], "Replace Si.cif?")
        @test AtomBond.save_from_browser!(b) && v.path == joinpath(dir, "Si.cif")
        @test save!(v) === true                            # now ours: saved in place
        b.format.i_selected[] = findfirst(f -> f[2] == ".xyz", AtomBond.SAVE_FORMATS)
        @test AtomBond.shown_text(b.name) == "Si.xyz"
        @test AtomBond.save_from_browser!(b) && formula(load_structure(joinpath(dir, "Si.xyz"))) == "Si8"
        @test AtomBond.with_format_extension("POSCAR", ".vasp") == "POSCAR"
        @test AtomBond.with_format_extension("Si.v2", ".cif") == "Si.v2.cif"
        @test AtomBond.save_extension("x.geom") == ".cell"
        b.format.i_selected[] = findfirst(f -> f[2] == ".png", AtomBond.SAVE_FORMATS)
        @test AtomBond.save_from_browser!(b) && isfile(joinpath(dir, "Si.png")) && v.path == joinpath(dir, "Si.xyz")
        # New structure: a new tab, unsaved, Save goes to Save as.
        d = AtomBond.new_structure_dialog(v)
        d.proto.i_selected[] = findfirst(==("Rocksalt"), AtomBond.prototype_names())
        @test AtomBond.shown_text(d.elements) == "Na Cl"
        AtomBond.settext!(d.elements, "Mg O"); AtomBond.settext!(d.abc, "4.21")
        d.create.clicks[] = 1
        @test formula(v.crystal) == "Mg4O4" && v.path === nothing && v.document[] == "Untitled — not saved"
        @test save!(v) isa AtomBond.FileBrowser
        n = ntabs(v)
        @test close_tab!(v; force=true) && ntabs(v) == n - 1
        @test v.path == joinpath(dir, "Si.xyz") && v.document[] == "Si.xyz"
        d.proto.i_selected[] = 1                                    # empty cell
        AtomBond.settext!(d.abc, "4 4 6"); AtomBond.settext!(d.angles, "90 90 120")
        d.create.clicks[] = 1
        @test natoms(v.crystal) == 0 && startswith(v.header[], "No atoms")
        # visualize() with no structure: an empty cell to build in.
        e = visualize(; show=false, panel=false)
        @test natoms(e.crystal) == 0 && e.path === nothing && lattice_parameters(e.crystal.lattice)[1:3] == (5.0, 5.0, 5.0)
        AtomBond.settext!(d.elements, "Si")
        d.create.clicks[] = 1
        @test startswith(v.status[], "Enter no elements")
        edit!(c -> add_atom(c, :Si, [0, 0, 0]), v, "add Si")
        @test formula(v.crystal) == "Si1"
        @test select_all!(v) !== nothing && length(v.selection) == 1
    end

    @testset "supercell dialog" begin
        v = visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false)
        @test "Supercell…" in [it.label for it in AtomBond.menubar_items(v, "Build")]
        d = AtomBond.supercell_dialog(v)
        txt(t) = AtomBond.shown_text(t)
        @test txt.(d.vecs) == ["2 0 0", "0 2 0", "0 0 2"] && startswith(d.preview[], "8× the cell: 64 atoms")
        d.repeats.displayed_string[] = "3x1x2"
        @test txt.(d.vecs) == ["3 0 0", "0 1 0", "0 0 2"]
        d.create.clicks[] = 1
        @test natoms(v.crystal) == 48 && startswith(v.status[], "Make 3×1×2 supercell")
        undo!(v)
        @test natoms(v.crystal) == 8 && startswith(d.preview[], "6× the cell: 48 atoms")
        d.repeats.displayed_string[] = "1 1 1"
        d.preset.i_selected[] = 2                                   # √2×√2 R45°
        @test txt.(d.vecs) == ["1 1 0", "-1 1 0", "0 0 1"] && occursin("a b c = 7.681 7.681 5.431", d.preview[])
        d.create.clicks[] = 1
        @test natoms(v.crystal) == 16
        undo!(v)
        AtomBond.settext!(d.vecs[2], "1 1 0")                        # dependent vectors
        @test startswith(d.preview[], "The new vectors must be")
        d.create.clicks[] = 1
        @test natoms(v.crystal) == 8 && startswith(v.status[], "Could not")
        AtomBond.settext!(d.vecs[2], "1 1")
        @test startswith(d.preview[], "Enter a′ b′ c′")
        # √3×√3 R30° of a hexagonal cell
        w = visualize(joinpath(STRUCTS, "GaN_wurtzite.vasp"); show=false)
        dw = AtomBond.supercell_dialog(w)
        dw.repeats.displayed_string[] = "1 1 1"; dw.preset.i_selected[] = 3
        dw.create.clicks[] = 1
        a, b, c, α, β, γ = lattice_parameters(w.crystal.lattice)
        @test natoms(w.crystal) == 12 && a ≈ b ≈ sqrt(3) * lattice_parameters(w.original.lattice)[1] && γ ≈ 120
    end

    @testset "shift and wrap atoms" begin
        g = read_poscar(joinpath(STRUCTS, "GaAs_conv.vasp"))
        a = lattice_parameters(g.lattice)[1]
        t = translate_atoms(g, [0.3, 0, 0]; fractional=true)
        @test t.frac[1] ≈ g.frac[1] + [0.3, 0, 0]
        @test translate_atoms(g, [0.3a, 0, 0]).frac ≈ t.frac                       # Cartesian Å
        @test translate_atoms(g, [0, 0, 0.5]; fractional=true, atoms=[2]).frac[1] == g.frac[1]
        @test atoms_outside_cell(t) == findall(f -> f[1] > 1, t.frac) && length(atoms_outside_cell(t)) == 2
        w = wrap(t)
        @test isempty(atoms_outside_cell(w)) && all(f -> all(0 .<= f .< 1), w.frac)
        @test wrap(t; atoms=[atoms_outside_cell(t)[1]]) |> atoms_outside_cell |> length == 1
        @test wrap(Crystal(cubic(4.0), [:Si], [[-1e-17, 1 - 1e-12, 0.5]])).frac[1] == [0, 0, 0.5]   # no 0.9999…
        @test atoms_outside_cell(Crystal(cubic(4.0), [:Si], [[1.0005, -0.0005, 0]])) == Int[]    # on a face

        # Drawn where they are; wrap_atoms=true draws them in the cell
        v = visualize(joinpath(STRUCTS, "GaAs_conv.vasp"); show=false)
        @test !any(i -> i.label == "Wrap atoms into cell" && i.kind == :item, AtomBond.menubar_items(v, "Build"))
        @test AtomBond.shift_atoms!(v, [0.3, 0, 0]; fractional=true)
        @test occursin("2 atoms are outside the cell", v.status[])
        xmax = maximum(at.position[1] for at in v.display.atoms)
        @test xmax ≈ 1.05a
        @test maximum(at.position[1] for at in AtomBond.build_display(v.crystal, ViewSettings(wrap_atoms=true)).atoms) <= a + 1e-6
        @test any(i -> i.label == "Wrap atoms into cell" && i.kind == :item && i.shortcut == "2 outside",
                  AtomBond.menubar_items(v, "Build"))
        # a whole lattice vector: same picture, atoms on faces stay on both faces
        si = visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false)
        n0 = length(si.display.atoms)
        AtomBond.shift_atoms!(si, [0, 0, 1]; fractional=true)
        @test count(at -> at.position[3] > 5.431 + 1e-6, si.display.atoms) == 8   # the 6 atoms off the z faces (8 copies)
        @test AtomBond.wrap_atoms!(si) && length(si.display.atoms) == n0
        @test !AtomBond.wrap_atoms!(si) && startswith(si.status[], "Nothing to wrap")

        # Dialog
        undo!(v)
        d = AtomBond.shift_dialog(v)
        txt(x) = AtomBond.shown_text(x)
        @test txt(d.vector) == "0 0 0" && d.atoms.i_selected[] == 1
        AtomBond.settext!(d.vector, "0.3 0 0")
        @test startswith(d.preview[], "Moves 8 atoms by") && endswith(d.preview[], "2 atoms will be outside the cell")
        d.coords.i_selected[] = 2                                    # → Cartesian, same shift
        @test isapprox(parse(Float64, split(txt(d.vector))[1]), 0.3a; atol=1e-4)
        d.coords.i_selected[] = 1
        @test txt(d.vector) == "0.3 0 0"
        d.wrap.active[] = true
        @test endswith(d.preview[], "(then wrapped into the cell)")
        d.apply.clicks[] = 1
        @test isempty(atoms_outside_cell(v.crystal)) && v.crystal.frac[1] ≈ [0.3, 0, 0]
        undo!(v)
        d.wrap.active[] = false
        # selected atoms only, the selection follows the moved atom
        v.selection = [(5, AtomBond.Shift(0, 0, 0))]; AtomBond.update_selection!(v)
        d.atoms.i_selected[] = 2
        AtomBond.settext!(d.vector, "0 0 1")
        @test startswith(d.preview[], "Moves 1 atom by")
        d.apply.clicks[] = 1
        @test v.crystal.frac[5] ≈ g.frac[5] + [0, 0, 1] && v.crystal.frac[1] == g.frac[1]
        @test AtomBond.selected_indices(v) == [5] && d.outside[] == "Selected atoms: 1 atom outside the cell"
        @test any(i -> i.label == "Wrap 1 atom into cell", AtomBond.menu_items(v))
        d.wrap_button.clicks[] = 1
        @test v.crystal.frac[5] ≈ g.frac[5]
        # selected atom → origin
        d.origin.clicks[] = 1
        @test d.atoms.i_selected[] == 1 && AtomBond.parse3(Float64, txt(d.vector)) ≈ -g.frac[5]
        d.apply.clicks[] = 1
        @test norm(v.crystal.frac[5]) < 1e-12 && spacegroup(v.crystal)[1] == 216
        AtomBond.settext!(d.vector, "0 0 0")
        d.apply.clicks[] = 1
        @test v.status[] == "Enter a non-zero shift"
    end

    @testset "rotate atoms" begin
        @test rotation_matrix([0, 0, 1], 90) * [1, 0, 0] ≈ [0, 1, 0]               # counterclockwise
        @test rotation_matrix([1, 1, 1], 120) * [1, 0, 0] ≈ [0, 1, 0]
        @test_throws ArgumentError rotation_matrix([0, 0, 0], 90)
        g = read_poscar(joinpath(STRUCTS, "GaAs_conv.vasp"))
        a = lattice_parameters(g.lattice)[1]
        # 90° about z through Ga1 puts As on the other tetrahedral sites
        S4 = wrap(rotate_atoms(g, 90, [0, 0, 1]; center=[0, 0, 0]))
        @test any(f -> f ≈ [0.75, 0.25, 0.25], S4.frac) && !any(f -> f ≈ [0.75, 0.25, 0.25], g.frac)
        @test spacegroup(rotate_atoms(g, 180, [0, 0, 1]; center=[0, 0, 0]))[1] == 216      # 2-fold axis
        r = rotate_atoms(g, 90, [0, 0, 1]; atoms=[5], center=[0, 0, 0])
        @test r.frac[5] ≈ [-0.25, 0.25, 0.25] && r.frac[1:4] == g.frac[1:4]
        # an image copy: the rotated copy is stored
        @test rotate_atoms(g, 0.0, [0, 0, 1]; atoms=[5], images=[(1, 0, 0)]).frac[5] ≈ g.frac[5] + [1, 0, 0]
        @test_throws ArgumentError rotate_atoms(g, 90, [0, 0, 1]; atoms=[5, 5])
        @test closest_contact(r, [5])[3] ≈ sqrt(3) / 4 * a                        # still bonded to Ga1
        hit = rotate_atoms(g, 180, [0, 0, 1]; atoms=[5], center=[0.375a, 0.375a, 0])   # → (½ ½ ¼)
        @test closest_contact(hit, [5])[3] ≈ a / 4
        @test closest_contact(g, [1]) |> last ≈ sqrt(3) / 4 * a

        # Water: rigid rotation about its centroid
        w = visualize(joinpath(STRUCTS, "water.xyz"); show=false)
        d = AtomBond.rotate_dialog(w)
        txt(x) = AtomBond.shown_text(x)
        @test (txt(d.axis), txt(d.angle), d.atoms.i_selected[]) == ("0 0 1", "90", 1)
        @test startswith(d.preview[], "Turns 3 atoms by 90° about [0 0 1]") && !occursin("cell does not turn", d.preview[])
        p0 = cartesian_positions(w.crystal)
        AtomBond.settext!(d.angle, "37")
        d.apply.clicks[] = 1
        p1 = cartesian_positions(w.crystal)
        @test startswith(w.status[], "Rotate 3 atoms by 37° about [0 0 1]")
        @test sum(p1) ≈ sum(p0) && all(norm(p1[i] - p1[j]) ≈ norm(p0[i] - p0[j]) for i in 1:3 for j in 1:3)
        @test all(p1[i][3] ≈ p0[i][3] for i in 1:3) && !(p1 ≈ p0)
        undo!(w)
        @test cartesian_positions(w.crystal) ≈ p0

        v = visualize(joinpath(STRUCTS, "GaAs_conv.vasp"); show=false)
        @test "Rotate atoms…" in [it.label for it in AtomBond.menubar_items(v, "Build")]
        d = AtomBond.rotate_dialog(v)
        @test occursin("The cell does not turn", d.preview[])
        # lattice axis [1 1 1], 120° through Ga1: a symmetry operation of zincblende
        d.axis_menu.i_selected[] = 2; AtomBond.settext!(d.axis, "1 1 1"); AtomBond.settext!(d.angle, "120")
        d.centre_menu.i_selected[] = 3; AtomBond.settext!(d.centre, "0 0 0")
        d.apply.clicks[] = 1
        @test spacegroup(v.crystal)[1] == 216 && occursin("[1 1 1] of the cell", v.status[])
        @test all(any(s == t && isapprox(f, h; atol=1e-9) for (t, h) in zip(g.species, g.frac))
                  for (s, f) in zip(v.crystal.species, wrap(v.crystal).frac))
        undo!(v)
        # selected atoms, axis through the first two clicked: atoms on the axis stay put
        v.selection = [(1, AtomBond.Shift(0, 0, 0)), (5, AtomBond.Shift(0, 0, 0)), (6, AtomBond.Shift(0, 0, 0))]
        AtomBond.update_selection!(v)
        d.atoms.i_selected[] = 2; d.axis_menu.i_selected[] = 3; AtomBond.settext!(d.angle, "120")
        @test startswith(d.preview[], "Turns 3 atoms by 120° about the Ga1–As5 axis")
        d.apply.clicks[] = 1
        @test v.crystal.frac[1] ≈ g.frac[1] && v.crystal.frac[5] ≈ g.frac[5] && !(v.crystal.frac[6] ≈ g.frac[6])
        # 120° about a Ga–As bond keeps As6 bonded to Ga1 (it moves to another neighbour site)
        @test isapprox(norm(cartesian_positions(v.crystal)[6]), sqrt(3) / 4 * a; atol=1e-9) ||
              isapprox(AtomBond.min_image_distance(v.crystal, 1, 6), sqrt(3) / 4 * a; atol=1e-9)
        @test AtomBond.selected_indices(v) == [1, 5, 6]
        v.selection = [(1, AtomBond.Shift(0, 0, 0))]; AtomBond.update_selection!(v)
        @test d.preview[] == "Select two atoms for the axis (the first two you clicked)"
        d.apply.clicks[] = 1
        @test startswith(v.status[], "Select two atoms")
        d.axis_menu.i_selected[] = 4                                   # viewing direction
        @test occursin("about the viewing direction", d.preview[])
        AtomBond.settext!(d.angle, "360")
        d.apply.clicks[] = 1
        @test startswith(v.status[], "Enter an angle")
    end

    @testset "exit" begin
        v = visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false)
        @test "Exit" in [it.label for it in AtomBond.menubar_items(v, "File")]
        @test !AtomBond.is_modified(v) && exit!(v) && v.status[] == "Closed"
        v = visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false)
        edit!(c -> AtomBond.delete_atoms(c, [1]), v, "delete 1 atom")
        @test !exit!(v) && startswith(v.status[], "Unsaved changes in Si_conv.vasp")
        edit!(c -> AtomBond.delete_atoms(c, [1]), v, "delete 1 atom")      # another edit: warn again
        @test !exit!(v)
        @test exit!(v)                                                    # second Exit closes
        w = visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false)
        edit!(c -> AtomBond.delete_atoms(c, [1]), w, "delete 1 atom"); undo!(w)
        @test exit!(w)                                                    # undone back to the saved state
        @test exit!(visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false); force=true)
    end

    @testset "structure library" begin
        lib = structure_library()
        groups = ["Semiconductors: group IV" => 3, "Semiconductors: III-V" => 12, "Oxides" => 10,
                  "Carbon: graphite and graphene" => 2, "Carbon: nanotubes" => 5, "Organic: hydrocarbons" => 11,
                  "Organic: halocarbons" => 12, "Organic: alcohols and ethers" => 3, "Organic: carbonyl compounds" => 5,
                  "Organic: nitrogen and sulfur compounds" => 10, "Organic: heterocycles" => 3, "Organic: radicals" => 4,
                  "Organometallic: metal alkyls" => 6, "Organometallic: metal carbonyls" => 3, "Inorganic: elements" => 9,
                  "Inorganic: hydrides" => 14, "Inorganic: oxides and sulfides" => 16, "Inorganic: oxoacids" => 4,
                  "Inorganic: nonmetal halides" => 26, "Inorganic: metal halides" => 11, "Inorganic: oxyhalides" => 6,
                  "Inorganic: radicals" => 6, "Carbon: fullerenes" => 1]
        @test length(lib) == 182 && [e.category for e in lib] == reduce(vcat, [fill(g, n) for (g, n) in groups])
        @test [e.name for e in lib[[1, 2, 3, 33, end]]] == ["Silicon", "Diamond", "Silicon-germanium Si0.5Ge0.5",
                                                             "Acetylene", "Buckminsterfullerene"]
        mols = structure_library(; kind=:molecule)                    # sorted by name within each group
        @test all(issorted(lowercase.([e.name for e in mols if e.category == g])) for (g, _) in groups[6:end])
        @test library_entry("SF6").use == "etch gas" && library_entry("Diborane").use == "dopant source" &&
              occursin("deposition precursor", library_entry("WF6").use) && library_entry("Benzene").use == ""
        # every molecule: the atoms of its formula, all connected by bonds, point group given
        count_atoms(f) = sum((isempty(m[1]) ? 1 : parse(Int, m[1]) for m in eachmatch(r"[A-Z][a-z]?(\d*)", f)); init=0)
        natoms_of(f) = count_atoms(replace(f, r"\(([^()]*)\)(\d*)" =>                 # Al(CH3)3 = Al C3 H9
                                           s -> (m = match(r"\(([^()]*)\)(\d*)", s); repeat(m[1], parse(Int, m[2])))))
        for e in structure_library(; kind=:molecule)
            c = load_library(e)
            n = natoms(c)
            seen, todo = Set(1), [1]                                  # flood fill along bonds
            # N2O4's N–N bond (1.782 Å) is longer than the default N–N limit; it needs a cutoff
            bonds = find_bonds(c; cutoffs=e.formula == "N2O4" ? Dict((:N, :N) => 1.85) : Dict{Tuple{Symbol,Symbol},Float64}())
            e.formula == "N2O4" && @test !any(b -> c.species[b.i] == c.species[b.j] == :N, find_bonds(c))
            while !isempty(todo)
                i = pop!(todo)
                for b in bonds, (x, y) in ((b.i, b.j), (b.j, b.i))
                    x == i && y ∉ seen && (push!(seen, y); push!(todo, y))
                end
            end
            @test (e.name, n, length(seen)) == (e.name, natoms_of(e.formula), n)
            @test occursin("point group", e.description)
        end
        @test occursin("triplet", library_entry("O2").note) && occursin("D∞h", library_entry("N2").description)
        # every crystal loads with the right composition and no overlapping atoms
        for e in structure_library(; kind=:crystal)
            c = load_library(e)
            @test (e.name, replace(formula(c), r"\d+" => "")) == (e.name, replace(e.formula, r"\d+" => ""))
            @test minimum(AtomBond.min_image_distance(c, i, j) for i in 1:natoms(c) for j in i+1:natoms(c)) > 1.35
        end
        sg(name) = spacegroup(load_library(name))[1]
        @test sg("GaAs") == 216 && sg("AlN") == 186 && sg("Al2O3") == 167 && sg("HfO2") == 14 && sg("SiO2") == 154 &&
              sg("SrTiO3") == 221 && sg("Ga2O3") == 12 && sg("Anatase (TiO2)") == 141 && sg("Rutile (TiO2)") == 136
        @test lattice_parameters(load_library("GaAs").lattice)[1] == 5.65325
        @test lattice_parameters(load_library("InSb").lattice)[1] == 6.4794
        @test all(e -> isfile(e.file) && !isempty(e.source) && !isempty(e.reference), lib)
        @test length(structure_library(; kind=:molecule)) == 150
        # crystals: published lattice constants and symmetry
        si = load_library("Silicon")
        @test si.title == "Silicon" && formula(si) == "Si8" && spacegroup(si)[1] == 227
        @test lattice_parameters(si.lattice)[1] == 5.430941
        @test lattice_parameters(load_library("Diamond").lattice)[1] == 3.566986
        sige = load_library("SiGe")                                   # by formula
        @test formula(sige) == "Si4Ge4" && spacegroup(sige)[1] == 216 && lattice_parameters(sige.lattice)[1] == 5.5373
        gan = load_library("GaN")
        @test formula(gan) == "Ga2N2" && spacegroup(gan)[1] == 186
        @test lattice_parameters(gan.lattice)[[1, 3]] == (3.189, 5.185)                    # bulk (Vurgaftman 2001)
        @test minimum(AtomBond.min_image_distance(gan, i, j) for i in 1:2 for j in 3:4) ≈ sqrt(3.189^2 / 3 + (0.5 - 0.3789)^2 * 5.185^2) atol=1e-6
        @test occursin("Paszkowicz", library_entry("GaN").reference)
        # molecules: experimental geometries from the literature, boxed
        bond(c, i, j) = norm(cartesian_positions(c)[i] - cartesian_positions(c)[j])
        w = load_library("H2O")
        @test w.pbc == (false, false, false) && isapprox(bond(w, 1, 2), 0.958; atol=1e-3)
        @test isapprox(AtomBond.bond_angle(cartesian_positions(w)[[2, 1, 3]]...), 104.48; atol=0.05)
        @test isapprox(bond(load_library("Methane"), 1, 2), 1.087; atol=1e-3)
        @test isapprox(bond(load_library("NH3"), 1, 2), 1.012; atol=1e-3)
        m = load_library("Methanol")
        @test formula(m) == "C1O1H4" && isapprox(bond(m, 1, 2), 1.427; atol=1e-3)
        # etch and dopant gases: shapes from the builders
        sf6 = load_library("SF6")
        @test natoms(sf6) == 7 && all(isapprox(bond(sf6, 1, k), 1.556; atol=1e-3) for k in 2:7)
        bf3 = cartesian_positions(load_library("BF3"))
        @test all(isapprox(AtomBond.bond_angle(bf3[i], bf3[1], bf3[j]), 120; atol=0.01) for (i, j) in ((2, 3), (3, 4), (2, 4)))
        pf5 = load_library("PF5")
        @test sort(round.([bond(pf5, 1, k) for k in 2:6]; digits=3)) == [1.532, 1.532, 1.532, 1.58, 1.58]
        clf3 = cartesian_positions(load_library("ClF3"))
        @test isapprox(AtomBond.bond_angle(clf3[2], clf3[1], clf3[3]), 87.483; atol=0.01) &&
              maximum(abs(p[3] - clf3[1][3]) for p in clf3) < 1e-6            # planar T shape
        @test occursin("r_g", library_entry("BCl3").note) && occursin("CCSD(T)", library_entry("CH3F").note)
        # precursors: metal alkyls and carbonyls
        tma2 = load_library("Al2(CH3)6")
        @test formula(tma2) == "Al2C6H18" && isapprox(bond(tma2, 1, 2), 2.619; atol=1e-3) &&
              isapprox(bond(tma2, 1, 3), 2.140; atol=1e-3)
        w6 = load_library("W(CO)6")
        @test natoms(w6) == 13 && all(isapprox(bond(w6, 1, k), 2.0586; atol=1e-3) for k in 2:7)
        @test isapprox(bond(load_library("Trimethylindium"), 1, 2), 2.161; atol=1e-3)
        al2cl6 = load_library("Al2Cl6")
        @test isapprox(bond(al2cl6, 1, 3), 2.250; atol=1e-3) && isapprox(bond(al2cl6, 1, 5), 2.061; atol=1e-3) &&
              isapprox(bond(al2cl6, 1, 2), 2.250 * sqrt(2); atol=1e-3)            # square Al2Cl2 ring
        sihcl3 = load_library("SiHCl3")
        @test isapprox(bond(sihcl3, 1, 2), 2.0118; atol=1e-3) && isapprox(bond(sihcl3, 1, 5), 1.4655; atol=1e-3)
        @test isapprox(bond(load_library("TiCl4"), 1, 2), 2.185; atol=1e-3) && occursin("Kimura", library_entry("TiCl4").reference)
        # organics: rings built from bonds and angles reproduce the papers' other values
        py = cartesian_positions(load_library("Pyridine"))
        @test isapprox(norm(py[1] - py[4]), 2.7980; atol=2e-4)                   # N1···C4, Table 4
        fu = cartesian_positions(load_library("Furan"))
        @test isapprox(norm(fu[3] - fu[4]), 1.4331; atol=2e-4) &&                # C3–C4
              isapprox(AtomBond.bond_angle(fu[2], fu[3], fu[4]), 105.992; atol=0.01)
        @test isapprox(norm(cartesian_positions(load_library("Thiophene"))[3] - cartesian_positions(load_library("Thiophene"))[4]),
                       1.4224; atol=2e-4)
        cb = cartesian_positions(load_library("Cyclobutane"))                   # puckered ring, 27.9°
        n1, n2 = cross(cb[4] - cb[1], cb[2] - cb[1]), cross(cb[4] - cb[3], cb[2] - cb[3])
        @test isapprox(180 - acosd(dot(n1, n2) / norm(n1) / norm(n2)), 27.9; atol=0.05)
        # carbon nanostructures: every atom has three C–C bonds near graphite's 1.4226 Å
        for name in ("Graphite", "Graphene", "Carbon nanotube (5,5)", "Carbon nanotube (6,5)", "C60")
            c = load_library(name)
            nb = zeros(Int, natoms(c))
            for b in find_bonds(c)
                nb[b.i] += 1; nb[b.j] += 1
                @test 1.40 < b.length < 1.46
            end
            @test all(==(3), nb)
        end
        @test load_library("Graphene").pbc == (true, true, false) && load_library("Carbon nanotube (6,5)").pbc == (false, false, true)
        @test natoms(load_library("Carbon nanotube (6,5)")) == 364
        # nanotube(n, m) for any indices
        for (n, m) in ((5, 5), (9, 0), (6, 5), (12, 3), (1, 1), (7, 0))
            c = nanotube(n, m)
            nb = zeros(Int, natoms(c))
            for b in find_bonds(c)
                nb[b.i] += 1; nb[b.j] += 1
            end
            @test (n, m, natoms(c)) == (n, m, 4 * (n^2 + n * m + m^2) ÷ gcd(2m + n, 2n + m))
            @test (n, m, all(==(3), nb)) == (n, m, true) || n + m < 4             # very thin tubes: odd contacts
        end
        t3 = nanotube(5, 5; periods=3, vacuum=12)
        @test natoms(t3) == 60 && t3.pbc == (false, false, true) && isapprox(lattice_parameters(t3.lattice)[3], 3 * 2.464; atol=1e-3)
        @test isapprox(lattice_parameters(t3.lattice)[1] - 12, 2 * 5 * 2.464 * sqrt(3) / 2π; atol=1e-3)   # box = diameter + vacuum
        bn = nanotube(5, 0; species=(:B, :N), bond=1.45)
        @test formula(bn) == "B10N10" && all(b -> bn.species[b.i] != bn.species[b.j], find_bonds(bn))
        @test cartesian_positions(load_library("Carbon nanotube (6,5)")) ≈ cartesian_positions(nanotube(6, 5)) atol=1e-5
        @test_throws ArgumentError nanotube(0, 0)
        @test_throws ArgumentError nanotube(4, 4; species=(:C, :Xx))
        @test AtomBond.nanotube_kind(9, 0) == ("zigzag", true) && AtomBond.nanotube_kind(6, 5) == ("chiral", false)
        dmso = load_library("Dimethyl sulfoxide")
        @test formula(dmso) == "S1O1C2H6" && isapprox(bond(dmso, 1, 3), 1.799; atol=1e-3)
        @test_throws ArgumentError load_library("Unobtainium")
        # your own library is added after the built-in one
        withenv("ATOMBOND_HOME" => mktempdir()) do
            dir = mkpath(AtomBond.user_library_dir())
            cp(joinpath(STRUCTS, "NaCl_conv.vasp"), joinpath(dir, "NaCl.vasp"))
            write(joinpath(dir, "library.toml"), "[[structure]]\nname = \"Rock salt\"\nkind = \"crystal\"\n" *
                                                 "formula = \"NaCl\"\nfile = \"NaCl.vasp\"\nsource = \"my file\"\n")
            @test structure_library()[end].name == "Rock salt" && formula(load_library("Rock salt")) == "Na4Cl4"
        end
        @test length(structure_library()) == 182

        # Window
        v = visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false)
        @test "Structure library…" in [it.label for it in AtomBond.menubar_items(v, "File")]
        @test "From library…" in [it.label for it in AtomBond.menubar_items(v, "Build")]
        @test "Nanotube…" in [it.label for it in AtomBond.menubar_items(v, "Build")]
        nt = AtomBond.nanotube_dialog(v)
        @test startswith(nt.preview[], "(6,5) chiral, semiconducting\ndiameter 7.4") && occursin("364 atoms", nt.preview[])
        AtomBond.settext!(nt.indices, "10 10"); AtomBond.settext!(nt.periods, "2")
        @test startswith(nt.preview[], "(10,10) armchair, metallic\n") && occursin("80 atoms in the cell (2 periods)", nt.preview[])
        AtomBond.settext!(nt.indices, "0 0")
        @test startswith(nt.preview[], "Enter the indices")
        AtomBond.settext!(nt.indices, "8 0"); AtomBond.settext!(nt.periods, "1"); AtomBond.settext!(nt.elements, "B N")
        tabs0 = ntabs(v)
        nt.create.clicks[] = 1
        @test ntabs(v) == tabs0 + 1 && formula(v.crystal) == "B16N16" && v.crystal.pbc == (false, false, true)
        close_tab!(v)
        d = AtomBond.library_dialog(v)
        @test startswith(d.info[], "Si, 8 atoms in the cell, Fd-3m (227)") && occursin("COD 9011998", d.info[])
        d.kind.i_selected[] = 3                                       # molecules
        @test length(d.entries[]) == 150 && startswith(d.info[], "C2H2, 4 atoms") && occursin("Group: Organic: hydrocarbons", d.info[])
        @test d.group.options[][1] == "All groups" && length(d.group.options[]) == 19 &&    # all + 18 molecule groups
              d.group.options[][2] == "Organic: hydrocarbons"
        d.group.i_selected[] = findfirst(==("Inorganic: hydrides"), d.group.options[])
        @test length(d.entries[]) == 14 && all(e -> e.category == "Inorganic: hydrides", d.entries[])
        d.group.i_selected[] = 1
        d.search.stored_string[] = "etch"                             # search: name, formula, group or use
        @test length(d.entries[]) == 15 && occursin("Used as: etch gas", d.info[])
        d.search.stored_string[] = "tungsten hexa"                    # all words must match
        @test [e.formula for e in d.entries[]] == ["W(CO)6", "WCl6", "WF6"]
        d.search.stored_string[] = "unobtainium"
        @test isempty(d.entries[]) && d.structure.options[] == ["(none)"] && d.info[] == "No structure matches"
        d.search.stored_string[] = "solvent"
        @test length(d.entries[]) == 8
        d.search.stored_string[] = ""
        @test length(d.entries[]) == 150
        d.structure.i_selected[] = findfirst(e -> e.name == "Methanol", d.entries[])
        @test occursin("Venkateswarlu", d.info[])
        d.open.clicks[] = 1
        @test formula(v.crystal) == "C1O1H4" && v.path === nothing && v.document[] == "Untitled — not saved"
        @test startswith(v.status[], "Opened Methanol from the library") && ntabs(v) == 2
        @test close_tab!(v; force=true)
        @test formula(v.crystal) == "Si8" && v.path !== nothing
        @test open_library!(v, "GaN") && formula(v.crystal) == "Ga2N2"
        # a molecule in a periodic box, for plane-wave codes
        d.boxed.active[] = true; AtomBond.settext!(d.vacuum, "12")
        d.structure.i_selected[] = findfirst(e -> e.name == "Water", d.entries[])
        d.open.clicks[] = 1
        @test formula(v.crystal) == "O1H2" && all(v.crystal.pbc) && occursin("12 Å vacuum box", v.status[])
        a, b, c = lattice_parameters(v.crystal.lattice)[1:3]
        @test a ≈ b ≈ c ≈ 12 + 2 * 0.7572 && kpoint_grid(v.crystal) == (1, 1, 1)
        @test length(v.plots) > 2                                     # the cell is drawn
    end

    @testset "constraints" begin
        g = load_library("GaN")
        @test !has_constraints(g) && g.cell_constraints == AtomBond.FREE_CELL && all(f -> f == (false, false, false), g.fixed)
        c = fix_atoms(g, [1])
        c = fix_atoms(c, [3], (false, false, true))
        c = with_constraints(c; cell=(0, 0, 1, 0, 0, 0))
        @test has_constraints(c) && c.fixed[[1, 3]] == [(true, true, true), (false, false, true)]
        @test describe_cell_constraints(c.cell_constraints) == "a, b fixed; c free; angles fixed"
        @test describe_cell_constraints((1, 1, 2, 0, 0, 0)) == "a = b together; c free; angles fixed"
        @test AtomBond.check_cell_constraints((5, 5, 9, 0, 0, 0)) == (1, 1, 2, 0, 0, 0)
        @test_throws ArgumentError AtomBond.check_cell_constraints((1, 2, 3, 1, 0, 0))   # length tied to angle
        # kept through edits that keep the atoms, repeated in supercells, cleared by free_atoms
        @test translate_atoms(c, [0.1, 0, 0]; fractional=true).fixed == c.fixed
        @test wrap(c).cell_constraints == (0, 0, 1, 0, 0, 0)
        @test AtomBond.delete_atoms(c, [2]).fixed == c.fixed[[1, 3, 4]]
        s = supercell(c, 2, 1, 1)
        @test count(all, s.fixed) == 2 && s.cell_constraints == (0, 0, 1, 0, 0, 0)
        @test supercell(c, [1 -1 0; 1 1 0; 0 0 1]).cell_constraints == AtomBond.FREE_CELL
        @test !any(any, free_atoms(c).fixed)
        # exporters
        @test AtomBond.cell_constraints_qe(c)[1:2] == ("vc-relax", "z")
        @test AtomBond.cell_constraints_qe(with_constraints(c; cell=(0, 0, 0, 0, 0, 0)))[1:2] == ("relax", nothing)
        @test AtomBond.cell_constraints_vasp(with_constraints(c; cell=(1, 1, 1, 0, 0, 0)))[1] == 7
        @test AtomBond.cell_constraints_vasp(c)[1] === nothing
        dir = mktempdir()
        for f in ("POSCAR", "x.cell", "x.in")
            r = load_structure(save_structure(joinpath(dir, f), c))
            @test (f, r.fixed) == (f, c.fixed)
        end
        @test load_structure(joinpath(dir, "x.cell")).cell_constraints == (0, 0, 1, 0, 0, 0)
        cell = read(joinpath(dir, "x.cell"), String)
        @test occursin("%BLOCK CELL_CONSTRAINTS\n  0 0 1\n  0 0 0", cell) && count("Ga     1", cell) == 3
        qe = read(joinpath(dir, "x.in"), String)
        @test occursin("calculation = 'vc-relax'", qe) && occursin("cell_dofree = 'z'", qe) && occursin("0 0 0\n", qe)
        @test !occursin("approximate", qe)                     # fixing z is exact here (c* along z)
        @test occursin("approximate", sprint(io -> write_qe(io, fix_atoms(g, [1], (true, false, false)))))
        @test occursin("Selective dynamics", read(joinpath(dir, "POSCAR"), String))
        @test !occursin("Selective", sprint(io -> write_poscar(io, g)))
        @test occursin("calculation = 'scf'", sprint(io -> write_qe(io, g)))
        paths = write_vasp_inputs(joinpath(dir, "vasp"), c)
        @test basename(paths[end]) == "INCAR.relax" && occursin("ISIF = 3   # not exact", read(paths[end], String))
        # a molecule in a periodic box
        w = load_library("H2O")
        b = periodic_box(w; vacuum=10, cubic=false)
        @test all(b.pbc)
        @test sort(collect(lattice_parameters(b.lattice)[1:3]) .- 10) ≈ [0, 0.5865, 2 * 0.7572] atol=1e-4
        @test cartesian_positions(b)[1] - cartesian_positions(b)[2] ≈ cartesian_positions(w)[1] - cartesian_positions(w)[2]
        @test_throws ArgumentError periodic_box(b)
        @test !any(as_molecule(b).pbc)

        # Window and actions (all undoable)
        v = visualize(joinpath(STRUCTS, "GaN_wurtzite.vasp"); show=false)
        @test "Cell and constraints…" in [it.label for it in AtomBond.menubar_items(v, "Build")]
        d = AtomBond.cell_dialog(v)
        @test startswith(d.info[], "Periodic crystal, 4 atoms") && occursin("density 6.", d.info[])
        @test d.preset.i_selected[] == 1 && occursin("QE: calculation = 'vc-relax', cell_dofree = 'all'", d.cell_preview[])
        d.preset.i_selected[] = 4                                     # only c relaxes
        @test AtomBond.shown_text(d.lengths) == "0 0 1" && occursin("cell_dofree = 'z'", d.cell_preview[])
        d.apply_cell.clicks[] = 1
        @test v.crystal.cell_constraints == (0, 0, 1, 0, 0, 0) && startswith(v.status[], "Set cell constraints")
        AtomBond.settext!(d.lengths, "1 1 2")
        @test d.preset.i_selected[] == 5
        AtomBond.settext!(d.angles, "1 0 0")
        d.apply_cell.clicks[] = 1
        @test startswith(v.status[], "Cannot use these cell constraints") && v.crystal.cell_constraints == (0, 0, 1, 0, 0, 0)
        undo!(v)
        @test v.crystal.cell_constraints == AtomBond.FREE_CELL && d.preset.i_selected[] == 1
        d.fix.clicks[] = 1
        @test startswith(v.status[], "Select atoms first")
        v.selection = [(1, AtomBond.Shift(0, 0, 0))]; AtomBond.update_selection!(v)
        d.axis_toggles[1].active[] = false; d.axis_toggles[2].active[] = false
        d.fix.clicks[] = 1
        @test v.crystal.fixed[1] == (false, false, true) && occursin("partly fixed: 1 atom", d.atoms_info[])
        @test endswith(v.info[], "fixed z")
        @test any(p -> p isa AtomBond.Makie.Scatter, v.plots)          # the + marker
        @test "Free 1 fixed atom" in [it.label for it in AtomBond.menu_items(v)]
        d.free_all.clicks[] = 1
        @test !has_constraints(v.crystal)
        d.box.clicks[] = 1                                            # periodic crystal → molecule
        @test !any(v.crystal.pbc) && startswith(d.info[], "Molecule")
        undo!(v)
        w = visualize(joinpath(STRUCTS, "water.xyz"); show=false)
        dw = AtomBond.cell_dialog(w)
        @test dw.box.label[] == "Put the molecule in a periodic box" && all(p -> p.visible[], w.triad.plots)  # axes shown
        # the molecule's box: one dashed grey line plot (no coloured a b c edges), and only with Unit cell on
        segs = filter(p -> p isa AtomBond.Makie.LineSegments, w.plots)
        col, box = segs[1].color[], w.settings.molecule_box_color
        @test length(segs) == 1 && (col.r, col.g, col.b) == (box.r, box.g, box.b) && segs[1].linestyle[] !== nothing
        w.settings.show_cell = false; redraw!(w)
        @test !any(p -> p isa AtomBond.Makie.LineSegments, w.plots)
        w.settings.show_cell = true; redraw!(w)
        dw.box.clicks[] = 1
        @test all(p -> p.visible[], w.triad.plots)
        @test all(w.crystal.pbc) && lattice_parameters(w.crystal.lattice)[1] ≈ lattice_parameters(w.crystal.lattice)[3]
    end

    @testset "bond length labels" begin
        v = visualize(load_library("Methanol"); show=false)
        texts(p) = p isa AtomBond.Makie.Text ? p[1][] isa AbstractVector ? p.text[] : String[] : String[]
        labels() = reduce(vcat, (texts(p) for p in [v.plots; v.overlay]); init=String[])
        @test isempty(labels())
        # C1–O2 1.427, O2–H6 0.956, C1–H3 1.096 among four selected atoms (H3–H6 are not bonded)
        v.selection = [(k, AtomBond.Shift(0, 0, 0)) for k in (1, 2, 3, 6)]; AtomBond.update_selection!(v)
        pts, t = AtomBond.bond_length_labels(v.display, v.settings; between=Set(v.selection))
        @test sort(t) == ["0.956", "1.096", "1.427"]
        @test issubset(t, labels())
        # two atoms: only the measured distance, no extra label
        v.selection = [(1, AtomBond.Shift(0, 0, 0)), (2, AtomBond.Shift(0, 0, 0))]; AtomBond.update_selection!(v)
        @test !("1.427" in labels()) && any(l -> startswith(l, "1.42"), labels())
        # all bonds (View menu), each once; then back
        @test "Label all atoms and bonds" in [it.label for it in AtomBond.menubar_items(v, "View")]
        set_bond_lengths!(v, :all)
        @test count(l -> occursin(r"^\d\.\d{3}$", l), labels()) == 5        # C–O, O–H, 3 × C–H
        @test "Label selected atoms only" in [it.label for it in AtomBond.menubar_items(v, "View")]
        set_bond_lengths!(v, :selected)
        @test count(l -> occursin(r"^\d\.\d{3}$", l), labels()) == 0
        # crystal: every bond and stub of zincblende GaAs is 2.448 Å; none when bonds are hidden
        g = visualize(joinpath(STRUCTS, "GaAs_conv.vasp"); show=false)
        set_bond_lengths!(g, :all)
        _, tg = AtomBond.bond_length_labels(g.display, g.settings)
        @test unique(tg) == ["2.448"] && length(tg) > 16
        g.settings.show_bonds = false; redraw!(g)
        @test !any(p -> p isa AtomBond.Makie.Text, g.plots)
        @test_throws ArgumentError set_bond_lengths!(g, :some)

        # Bond angles: off by default; toggle in the View menu; same scope as the lengths
        degrees(l) = occursin(r"^\d+\.\d°$", l)
        set_bond_lengths!(v, :selected)
        @test !v.settings.show_bond_angles && "Show bond angles" in [it.label for it in AtomBond.menubar_items(v, "View")]
        set_bond_angles!(v, true)
        @test "Hide bond angles" in [it.label for it in AtomBond.menubar_items(v, "View")]
        v.selection = [(k, AtomBond.Shift(0, 0, 0)) for k in (1, 2, 3, 6)]; AtomBond.update_selection!(v)
        # among C1 O2 H3 H6: O2–C1–H3 at C and C1–O2–H6 at O (4 selected show the dihedral, not an angle)
        @test sort(filter(degrees, labels())) == sort(AtomBond.bond_angle_labels(v.display, v.settings;
                                                                                  between=Set(v.selection))[2])
        @test count(degrees, labels()) == 2
        v.selection = [(k, AtomBond.Shift(0, 0, 0)) for k in (1, 2, 3)]; AtomBond.update_selection!(v)
        @test count(degrees, labels()) == 0                           # 3 selected: the measured angle only
        set_bond_lengths!(v, :all)
        _, ang = AtomBond.bond_angle_labels(v.display, v.settings)
        @test length(ang) == 7 && count(degrees, labels()) == 7        # 6 at C, 1 at O
        @test all(a -> 100 < parse(Float64, chop(a)) < 115, ang)
        set_bond_angles!(v, false)
        @test count(degrees, labels()) == 0
        # silicon: every angle tetrahedral; CO2: 180° placed beside C
        si = visualize(load_library("Silicon"); settings=ViewSettings(bond_lengths=:all, show_bond_angles=true), show=false)
        @test unique(AtomBond.bond_angle_labels(si.display, si.settings)[2]) == ["109.5°"]
        co2 = visualize(load_library("CO2"); settings=ViewSettings(bond_lengths=:all, show_bond_angles=true), show=false)
        p, t = AtomBond.bond_angle_labels(co2.display, co2.settings)
        @test t == ["180.0°"] && norm(AtomBond.Vec3(p[1]) - cartesian_positions(co2.crystal)[findfirst(==(:C), co2.crystal.species)]) > 0.3

        # Dihedral angles: toggle, values, placement, scope
        dih(w) = AtomBond.dihedral_labels(w.display, w.settings)
        e = visualize(load_library("Ethane"); settings=ViewSettings(bond_lengths=:all), show=false)
        @test !e.settings.show_dihedrals && "Show dihedral angles" in [it.label for it in AtomBond.menubar_items(e, "View")]
        set_dihedrals!(e, true)
        @test "Hide dihedral angles" in [it.label for it in AtomBond.menubar_items(e, "View")]
        p, t = dih(e)
        @test sort(t) == sort(["60.0°", "60.0°", "60.0°", "-60.0°", "-60.0°", "-60.0°", "180.0°", "180.0°", "180.0°"])
        @test length(unique(p)) == 9                                    # no two labels in the same place
        mid = sum(cartesian_positions(e.crystal)[1:2]) / 2
        @test all(q -> norm(AtomBond.Vec3(q) - mid) > 0.3, p)            # none on the C–C length label
        elabels() = reduce(vcat, (p.text[] for p in [e.plots; e.overlay] if p isa AtomBond.Makie.Text); init=String[])
        @test count(endswith("°"), elabels()) == 9
        @test only(dih(visualize(load_library("H2O2"); settings=ViewSettings(bond_lengths=:all, show_dihedrals=true),
                                 show=false))[2]) == "119.8°"
        @test isempty(dih(visualize(load_library("CO2"); settings=ViewSettings(bond_lengths=:all, show_dihedrals=true),
                                    show=false))[2])                     # straight: undefined
        si = visualize(load_library("Silicon"); settings=ViewSettings(bond_lengths=:all, show_dihedrals=true), show=false)
        @test sort(unique(dih(si)[2])) == ["-60.0°", "180.0°", "60.0°"]     # staggered; -180 is shown as 180
        # selected atoms: 5 or more (4 show the measured dihedral)
        set_bond_lengths!(e, :selected)
        H(k) = (k, AtomBond.Shift(0, 0, 0))
        e.selection = [H(3), H(1), H(2), H(6)]; AtomBond.update_selection!(e)
        @test count(endswith("°"), elabels()) == 1 && startswith(e.info[], "H3 – C1 – C2 – H6   dihedral")
        e.selection = [H(3), H(1), H(2), H(6), H(7)]; AtomBond.update_selection!(e)
        @test count(l -> occursin(r"^-?\d+\.\d°$", l), elabels()) == 2  # H3–C1–C2–H6 and H3–C1–C2–H7
        set_dihedrals!(e, false)
        @test count(l -> occursin(r"^-?\d+\.\d°$", l), elabels()) == 0

        # Atom coordinates: selected atoms (even one), or all; fractional for crystals, Cartesian for molecules
        coords(w) = filter(l -> startswith(l, "("), reduce(vcat, (p.text[] for p in [w.plots; w.overlay]
                                                                    if p isa AtomBond.Makie.Text); init=String[]))
        si = visualize(load_library("Silicon"); show=false)
        @test !si.settings.show_coordinates && "Show atom coordinates" in [it.label for it in AtomBond.menubar_items(si, "View")]
        set_coordinates!(si, true)
        @test isempty(coords(si))                                       # nothing selected
        si.selection = [(1, AtomBond.Shift(0, 0, 0))]; AtomBond.update_selection!(si)
        @test coords(si) == ["(0.000, 0.000, 0.000)"]
        @test "Coordinates: use Cartesian (Å)" in [it.label for it in AtomBond.menubar_items(si, "View")]
        set_coordinates!(si, true; system=:cartesian)
        si.selection = [(1, AtomBond.Shift(1, 0, 0))]; AtomBond.update_selection!(si)   # the copy on the far face
        @test coords(si) == ["(5.431, 0.000, 0.000)"] && "Coordinates: use fractional" in
              [it.label for it in AtomBond.menubar_items(si, "View")]
        set_coordinates!(si, true; system=:fractional)
        @test coords(si) == ["(1.000, 0.000, 0.000)"]
        set_bond_lengths!(si, :all)
        @test length(coords(si)) == length(si.display.atoms) == 18
        set_coordinates!(si, false)
        @test isempty(coords(si))
        @test_throws ArgumentError set_coordinates!(si, true; system=:polar)
        w = visualize(load_library("Water"); settings=ViewSettings(bond_lengths=:all, show_coordinates=true), show=false)
        @test AtomBond.coordinate_system(w.crystal, w.settings) === :cartesian
        @test sort(coords(w)) == sort([AtomBond.Printf.@sprintf("(%.3f, %.3f, %.3f)", r...) for r in cartesian_positions(w.crystal)])
    end

    @testset "atom charges" begin
        F = joinpath(STRUCTS, "formats", "GaAs_relax.castep")
        c = load_structure(F)
        @test structure_format(F) === :castep_out && has_charges(c) && c.charges == [0.46, -0.46]
        @test c.frac[2] ≈ [0.25, 0.25, 0.25] && lattice_parameters(c.lattice)[1] ≈ 2.826625 * sqrt(2)   # final step
        @test read_castep_output(F; charges=:hirshfeld).charges == [0.12, -0.12]
        @test !has_charges(read_castep_output(F; charges=:none))
        # Hirshfeld only (no Mulliken table) is used as the fallback
        txt = read(F, String)
        noM = replace(txt, r"Atomic Populations \(Mulliken\).*?(?=     Hirshfeld)"s => "")
        @test read_castep_output(IOBuffer(noM)).charges == [0.12, -0.12]
        # "Charge(e)" without a space (some CASTEP versions)
        @test read_castep_output(IOBuffer(replace(txt, "Charge (e)" => "Charge(e)"))).charges == [0.46, -0.46]
        # kept through edits; supercells repeat; a new element loses it; added atoms have none
        s = supercell(c, 2, 1, 1)
        @test s.charges == [0.46, -0.46, 0.46, -0.46] || sort(s.charges) == [-0.46, -0.46, 0.46, 0.46]
        @test isnan(substitute(c, [1], :Al).charges[1]) && substitute(c, [1], :Ga).charges[1] == 0.46
        @test translate_atoms(c, [0.1, 0, 0]; fractional=true).charges == c.charges
        @test AtomBond.delete_atoms(c, [1]).charges == [-0.46]
        @test isnan(add_atom(c, :H, [0.5, 0.5, 0.5]).charges[end])
        # files: extended XYZ and LAMMPS keep them; read from PDB formal charges
        dir = mktempdir()
        for f in ("x.xyz", "data.lmp")
            @test (f, load_structure(save_structure(joinpath(dir, f), c)).charges) == (f, [0.46, -0.46])
        end
        @test occursin("Atoms # charge", read(joinpath(dir, "data.lmp"), String))
        @test !occursin("charge", sprint(io -> write_xyz(io, load_library("Silicon"))))
        pdb = "CRYST1    4.000    4.000    4.000  90.00  90.00  90.00 P 1           1\n" *
              "HETATM    1 Na   MOL A   1       0.000   0.000   0.000  1.00  0.00          NA1+\n" *
              "HETATM    2 Cl   MOL A   1       2.000   2.000   2.000  1.00  0.00          CL1-\nEND\n"
        @test read_pdb(IOBuffer(pdb)).charges == [1.0, -1.0]
        @test AtomBond.charge_column("Lattice=\"1 0 0 0 1 0 0 0 1\" Properties=species:S:1:pos:R:3:forces:R:3:initial_charges:R:1") == 8

        # Display
        v = visualize(F; show=false)
        qlabels() = filter(l -> occursin(r"^[+−]?\d\.\d\d$", l), reduce(vcat, (p.text[] for p in [v.plots; v.overlay]
                                                                              if p isa AtomBond.Makie.Text); init=String[]))
        @test !v.settings.show_charges && "Show atom charges" in [it.label for it in AtomBond.menubar_items(v, "View")]
        set_charge_labels!(v, true)
        @test isempty(qlabels())
        v.selection = [(1, AtomBond.Shift(0, 0, 0))]; AtomBond.update_selection!(v)
        @test qlabels() == ["+0.46"] && endswith(v.info[], "charge +0.46 e")
        set_bond_lengths!(v, :all)
        @test sort(unique(qlabels())) == ["+0.46", "−0.46"] && length(qlabels()) == length(v.display.atoms)
        # bond lengths can be hidden on their own
        @test "Hide bond lengths" in [it.label for it in AtomBond.menubar_items(v, "View")]
        lens() = count(l -> occursin(r"^\d\.\d{3}$", l), reduce(vcat, (p.text[] for p in v.plots if p isa AtomBond.Makie.Text);
                                                                  init=String[]))
        @test lens() > 0
        set_bond_length_labels!(v, false)
        @test lens() == 0 && !isempty(qlabels())
        set_charge_labels!(v, false)
        @test isempty(qlabels())
        w = visualize(load_library("Water"); show=false)
        set_charge_labels!(w, true)
        @test startswith(w.status[], "Cannot show charges: this structure has none")
        @test AtomBond.format_charge(-0.001) == "0.00" && AtomBond.format_charge(1.234) == "+1.23"
    end

    @testset "magnetic moments" begin
        F = joinpath(STRUCTS, "formats", "Cr_afm.cell")
        c = load_structure(F)
        @test has_magmoms(c) && c.magmoms == [1.0, -1.0]
        # symmetry sees the magnetic order: AFM bcc is simple cubic with 2 atoms
        @test spacegroup(c)[1] == 221 && spacegroup(with_magmoms(c, nothing))[1] == 229
        @test natoms(primitive_cell(c)) == 2 && sort(primitive_cell(c).magmoms) == [-1.0, 1.0]
        @test natoms(primitive_cell(with_magmoms(c, nothing))) == 1
        # kept through edits and supercells; a new element loses it
        @test supercell(c, 2, 1, 1).magmoms == [1.0, -1.0, 1.0, -1.0] || sort(supercell(c, 2, 1, 1).magmoms) == [-1, -1, 1, 1]
        @test isnan(substitute(c, [1], :Fe).magmoms[1]) && AtomBond.delete_atoms(c, [1]).magmoms == [-1.0]
        @test wrap(translate_atoms(c, [0.5, 0, 0]; fractional=true)).magmoms == c.magmoms
        # files: CASTEP .cell SPIN=, extended XYZ; VASP gets ISPIN and MAGMOM in POSCAR order
        dir = mktempdir()
        for f in ("x.cell", "x.xyz")
            @test (f, load_structure(save_structure(joinpath(dir, f), c)).magmoms) == (f, [1.0, -1.0])
        end
        @test occursin("SPIN=-1.0000", read(joinpath(dir, "x.cell"), String))
        @test !occursin("SPIN", sprint(io -> write_castep(io, load_library("Silicon"))))
        paths = write_vasp_inputs(joinpath(dir, "vasp"), supercell(c, 2, 1, 1))
        @test basename(paths[end]) == "INCAR.spin" && occursin("ISPIN = 2\nMAGMOM = 1 -1 1 -1", read(paths[end], String))
        m = with_magmoms(load_library("GaAs"), [0, 0, 0, 0, 2, 2, 2, 0])
        incar = read(write_vasp_inputs(joinpath(dir, "v2"), m)[end], String)
        @test occursin("MAGMOM = 4*0 3*2 0", incar)
        # CASTEP output: Mulliken "Spin (hbar/2)" in both layouts (one line; up/dn lines)
        head = "     Atomic Populations (Mulliken)\n     -----------------------------\n"
        one = head * "Species          Ion     s       p       d       f      Total   Charge (e)   Spin (hbar/2)\n" *
              "========\n  Cr  1  2.20  6.50  4.90  0.00  13.60  0.40   1.10\n" *
              "  Cr  2  2.20  6.50  4.90  0.00  13.60  0.40  -1.10\n========\n"
        two = head * "Species          Ion Spin      s       p       d       f      Total   Charge(e)   Spin(hbar/2)\n" *
              "========\n  Cr  1  up:  1.10  3.25  3.00  0.00  7.35\n        dn:  1.10  3.25  1.90  0.00  6.25  0.40  1.10\n" *
              "  Cr  2  up:  1.10  3.25  1.90  0.00  6.25\n        dn:  1.10  3.25  3.00  0.00  7.35  0.40  -1.10\n========\n"
        for t in (one, two)
            L = split(t, '\n')
            @test AtomBond.castep_population(L, r"Atomic Populations \(Mulliken\)"i, "Spin") == [(:Cr, 1, 1.1), (:Cr, 2, -1.1)]
            @test AtomBond.castep_population(L, r"Atomic Populations \(Mulliken\)"i, "Charge") == [(:Cr, 1, 0.4), (:Cr, 2, 0.4)]
        end
        out = read(joinpath(STRUCTS, "formats", "GaAs_relax.castep"), String)
        @test !has_magmoms(read_castep_output(IOBuffer(out)))                          # not spin-polarized
        mag = replace(out, r"Species          Ion     s.*?(?=\n\n     Hirshfeld)"s =>
              "Species          Ion     s       p       d       f      Total   Charge (e)   Spin (hbar/2)\n" *
              "=====\n  Ga  1  1.30  1.24  10.00  0.00  12.54  0.46  0.20\n  As  1  1.66  3.80  0.00  0.00  5.46  -0.46  -0.20\n=====")
        r = read_castep_output(IOBuffer(mag))
        @test r.magmoms == [0.2, -0.2] && r.charges == [0.46, -0.46]

        # Display
        v = visualize(F; show=false)
        mlabels() = filter(l -> endswith(l, " μB"), reduce(vcat, (p.text[] for p in [v.plots; v.overlay]
                                                                    if p isa AtomBond.Makie.Text); init=String[]))
        @test !v.settings.show_magmoms && "Show magnetic moments" in [it.label for it in AtomBond.menubar_items(v, "View")]
        set_magmom_labels!(v, true)
        @test isempty(mlabels())
        v.selection = [(2, AtomBond.Shift(0, 0, 0))]; AtomBond.update_selection!(v)
        @test mlabels() == ["−1.00 μB"] && endswith(v.info[], "moment −1.00 μB")
        set_bond_lengths!(v, :all)
        @test sort(unique(mlabels())) == ["+1.00 μB", "−1.00 μB"] && length(mlabels()) == length(v.display.atoms)
        set_magmom_labels!(v, false)
        @test isempty(mlabels())
        w = visualize(load_library("Silicon"); show=false)
        set_magmom_labels!(w, true)
        @test startswith(w.status[], "Cannot show magnetic moments")
    end

    @testset "more file formats" begin
        F = joinpath(STRUCTS, "formats")
        sg(f) = spacegroup(load_structure(joinpath(F, f)))[1]
        expect = Dict("Si_scf.pwi" => ("Si2", 227), "GaN_relax.pwo" => ("Ga2N2", 186), "geometry.in" => ("Ga2N2", 186),
                      "Si.fdf" => ("Si2", 227), "GaAs.abi" => ("Ga1As1", 216), "vasprun.xml" => ("Na1Cl1", 225),
                      "ZnO.xsf" => ("Zn2O2", 186), "NaCl.pdb" => ("Na4Cl4", 225), "Si_hex.lmp" => ("Si2", 227),
                      "TiO2_rutile.car" => ("Ti2O4", 136), "Si.xsd" => ("Si2", 227), "GaN_airss.res" => ("Ga2N2", 186),
                      "Si_md.XDATCAR" => ("Si8", 35), "water.cube" => ("O1H2", 25),
                      "GaAs_relax.castep" => ("Ga1As1", 216), "Cr_afm.cell" => ("Cr2", 221))
        @test sort(collect(keys(expect))) == sort(readdir(F))
        for (f, (fml, n)) in expect
            c = load_structure(joinpath(F, f))
            @test (f, formula(c), spacegroup(c)[1]) == (f, fml, n)
        end
        # Lattice constants converted from Bohr / alat / scaled units
        @test isapprox(cell_volume(load_structure(joinpath(F, "Si_scf.pwi")).lattice), 5.431^3 / 4; rtol=1e-4)
        @test isapprox(cell_volume(load_structure(joinpath(F, "GaAs.abi")).lattice), 5.653^3 / 4; rtol=1e-4)
        @test cell_volume(load_structure(joinpath(F, "vasprun.xml")).lattice) ≈ 2 * 2.82^3     # finalpos, not initial
        # pw.x output: the final coordinates, not the first BFGS step
        q = load_structure(joinpath(F, "GaN_relax.pwo"))
        @test q.frac[3][3] ≈ 0.377
        @test load_structure(joinpath(F, "Si_md.XDATCAR")).frac[1][1] ≈ 0.01             # last frame
        @test read_xdatcar(joinpath(F, "Si_md.XDATCAR"); frame=1).frac[1][1] == 0
        w = load_structure(joinpath(F, "water.cube"))
        @test isapprox(minimum(norm(cartesian_positions(w)[1] - cartesian_positions(w)[k]) for k in 2:3), 0.957; atol=0.01)
        # QE: ibrav variants and units
        qe(sys, cards) = read_qe(IOBuffer("&control\n/\n&system\n$sys\n/\n&electrons\n/\n$cards"))
        c = qe("ibrav=1, A=4.0, nat=1, ntyp=1", "ATOMIC_SPECIES\nPo 209 Po.UPF\nATOMIC_POSITIONS crystal\nPo 0 0 0\n")
        @test cell_volume(c.lattice) ≈ 64
        c = qe("ibrav=4, celldm(1)=6.0263, celldm(3)=1.626, nat=2, ntyp=1",
               "ATOMIC_POSITIONS {crystal}\nMg 0.3333333 0.6666667 0.25\nMg 0.6666667 0.3333333 0.75\n")
        @test spacegroup(c)[1] == 194 && lattice_parameters(c.lattice)[6] ≈ 120
        c = qe("ibrav = 0\n nat = 1, ntyp = 1", "CELL_PARAMETERS angstrom\n 3 0 0\n 0 4 0\n 0 0 5\nATOMIC_POSITIONS angstrom\nFe1 1.5 2.0 2.5\n")
        @test c.species == [:Fe] && c.frac[1] ≈ [0.5, 0.5, 0.5]
        c = qe("ibrav=14, celldm(1)=10, celldm(2)=1.1, celldm(3)=1.2, celldm(4)=0.1, celldm(5)=0.2, celldm(6)=0.3, nat=1",
               "ATOMIC_POSITIONS crystal\nSi 0 0 0\n")
        p = lattice_parameters(c.lattice)
        @test p[1] ≈ 10 * AtomBond.BOHR && p[2] ≈ 11 * AtomBond.BOHR && collect(p[4:6]) ≈ acosd.([0.1, 0.2, 0.3])
        c = qe("ibrav=3, celldm(1)=5.42, nat=1", "ATOMIC_POSITIONS crystal\nFe 0 0 0\n")
        @test spacegroup(c)[1] == 229
        @test_throws ArgumentError qe("ibrav=13, celldm(1)=5, nat=1", "ATOMIC_POSITIONS crystal\nFe 0 0 0\n")
        @test_throws ArgumentError read_qe(IOBuffer("&control\n/\n"))
        # SHELX with centring + inversion generates the full cell
        r = read_res(IOBuffer("TITL NaCl\nCELL 0.71 5.64 5.64 5.64 90 90 90\nLATT 4\nSFAC Na Cl\n" *
                              "Na1 1 0 0 0 11.0 0.01\nCl1 2 0.5 0.5 0.5 11.0 0.01\nHKLF 4\nEND\n"))
        @test formula(r) == "Na4Cl4" && spacegroup(r)[1] == 225
        r = read_res(IOBuffer("TITL x\nCELL 1 4 4 4 90 90 90\nLATT -1\nSYMM -x,-y,z\nSFAC C\nC1 1 0.1 0.2 0.3 1.0\nEND\n"))
        @test natoms(r) == 2 && r.frac[2] ≈ [0.9, 0.8, 0.3]
        # Molecules (no cell) are boxed and non-periodic
        m = read_pdb(IOBuffer("ATOM      1  O   HOH A   1       0.000   0.000   0.117  1.00  0.00           O\n" *
                              "ATOM      2  H1  HOH A   1       0.000   0.757  -0.469  1.00  0.00           H\nEND\n"))
        @test formula(m) == "O1H1" && m.pbc == (false, false, false)
        m = read_aims(IOBuffer("atom 0 0 0 C\natom 0 0 1.09 H\n"))
        @test !any(m.pbc) && natoms(m) == 2
        m = read_xsf(IOBuffer("MOLECULE\nATOMS\n6 0 0 0\n1 0 0 1.09\n"))
        @test !any(m.pbc) && m.species == [:C, :H]
        # ABINIT: angdeg and xangst, n*x repetition
        a = read_abinit(IOBuffer("acell 3*3.0 angstrom\nangdeg 3*90\nnatom 2 ntypat 1 typat 2*1 znucl 26\n" *
                                 "xangst 0 0 0  1.5 1.5 1.5\n"))
        @test formula(a) == "Fe2" && spacegroup(a)[1] == 229
        # LAMMPS charge style, element from the mass
        l = read_lammps_data(IOBuffer("x\n\n1 atoms\n1 atom types\n0 3 xlo xhi\n0 3 ylo yhi\n0 3 zlo zhi\n\n" *
                                      "Masses\n\n1 15.999\n\nAtoms\n\n1 1 -2.0 1.5 1.5 1.5\n"))
        @test l.species == [:O] && l.frac[1] ≈ [0.5, 0.5, 0.5]

        # Writers: write, read back, same cell and atoms (low-symmetry triclinic cell)
        t = Crystal(Lattice(4.1, 5.2, 6.3, 80, 95, 110), [:Ga, :N, :O], [[0.1, 0.2, 0.3], [0.55, 0.6, 0.7], [0.9, 0.05, 0.45]];
                    title="test")
        dir = mktempdir()
        for (name, tol) in [("t.xsf", 1e-8), ("geometry.in", 1e-8), ("t.pdb", 2e-3), ("t.lmp", 1e-8), ("t.car", 1e-7),
                            ("t.res", 1e-6), ("t.cell", 1e-8), ("t.in", 1e-8), ("t.cif", 1e-5), ("t.xyz", 1e-6)]
            path = save_structure(joinpath(dir, name), t)
            u = load_structure(path)
            @test (name, u.species) == (name, t.species)
            @test all(isapprox.(lattice_parameters(u.lattice), lattice_parameters(t.lattice); atol=tol * 100))
            @test all(norm((f - g) .- round.(f - g)) < tol * 100 for (f, g) in zip(u.frac, t.frac))
        end
        # Non-periodic structures keep pbc in XSF and FHI-aims
        mol = read_xyz(IOBuffer("2\nCO\nC 0 0 0\nO 0 0 1.128\n"))
        for name in ("m.xsf", "geometry.in")
            u = load_structure(save_structure(joinpath(dir, name), mol))
            @test !any(u.pbc) && isapprox(norm(cartesian_positions(u)[2] - cartesian_positions(u)[1]), 1.128; atol=1e-6)
        end
        # Browser / Save as know the formats
        @test structure_format("run/XDATCAR") === :xdatcar && structure_format("vasprun.xml") === :vasprun
        @test structure_format("data.si") === :lammps && structure_format("geometry.in.next_step") === :aims
        @test all(f -> haskey(AtomBond.FORMAT_LABELS, f), unique(values(AtomBond.STRUCTURE_EXTENSIONS)))
        @test AtomBond.with_format_extension("GaN.vasp", "geometry.in") == "geometry.in"
        @test AtomBond.save_extension("x.xsd") == ".car" && AtomBond.save_extension("geometry.in") == "geometry.in"
    end

    @testset "click feedback" begin
        M = AtomBond.Makie
        v = visualize(joinpath(STRUCTS, "Si_conv.vasp"); show=false)
        AtomBond.GLMakie.activate!(ssao=false)
        frames(n=12) = for _ in 1:n; M.colorbuffer(v.figure); sleep(0.03); end   # drawn frames drive the effects
        frames(1)
        # A button flashes until it has been drawn and FLASH_TIME has passed.
        b = filter(x -> x isa M.Button && x.label[] == "Reset", v.figure.content)[1]
        normal = b.buttoncolor[]
        b.clicks[] += 1
        @test b.buttoncolor[] == AtomBond.FLASH_COLOR
        frames()
        @test b.buttoncolor[] == normal
        # A chosen menu item flashes; the menu closes at once but disappears after the flash.
        AtomBond.open_menubar_menu!(v, 4)
        m = v.dialogs[:menu]
        AtomBond.choose_menu_item!(v, 1)
        @test !AtomBond.menu_isopen(v) && m.scene.visible[] && m.flash[]
        frames()
        @test !m.scene.visible[]
        AtomBond.open_menubar_menu!(v, 4)
        @test !m.flash[] && m.scene.visible[]
        AtomBond.close_context_menu!(v)
        # The menu stays on top of plots added after it (redraws of the 3D view, a rebuilt tab bar):
        # GLMakie draws plots by z offset, so every menu plot must be above every other plot.
        redraw!(v); open_in_tab!(v, load_library("Water")); close_tab!(v)
        AtomBond.open_menubar_menu!(v, 4)
        zof(p) = AtomBond.Makie.transformationmatrix(p)[][3, 4]
        # plots that are drawn (visible, in visible scenes); Makie's hidden helpers may sit higher
        shown(s) = s.visible[] ? [filter(p -> p.visible[], s.plots); reduce(vcat, (shown(c) for c in s.children); init=Any[])] : Any[]
        others = filter(p -> !(p in m.scene.plots), shown(v.figure.scene))
        @test !isempty(m.scene.plots) && minimum(zof, m.scene.plots) > maximum(zof, others)
        AtomBond.close_context_menu!(v)
        # Ripple in the 3D view, cleared after RIPPLE_TIME.
        AtomBond.click_ripple!(v, (300.0, 300.0), :atom)
        r = v.dialogs[:ripple]
        @test !isnan(r.points[][1][1]) && r.color[].alpha > 0.5
        frames(25)
        @test isnan(r.points[][1][1])
        # Off switch
        @test click_feedback!(flash=false, sound=false) == (flash=false, sound=false)
        b.clicks[] += 1
        @test b.buttoncolor[] == normal
        AtomBond.click_ripple!(v, (300.0, 300.0))
        @test isnan(r.points[][1][1]) && !AtomBond.click_sound()
        click_feedback!(flash=true, sound=true)
        # The click sound file is a valid 25 ms WAV.
        wav = read(AtomBond.click_wav())
        @test String(wav[1:4]) == "RIFF" && String(wav[9:12]) == "WAVE" && length(wav) == 44 + 2 * round(Int, 0.025 * 44100)
    end

    @testset "less common file variants" begin
        # QE Bravais lattices: one atom per cell gives the lattice's own space group,
        # which checks the centring of each ibrav; the volume checks the lengths.
        qe(sys) = read_qe(IOBuffer("&system\nnat=1, ntyp=1, $sys\n/\nATOMIC_POSITIONS crystal\nFe 0 0 0\n"))
        a, b, c = 10AtomBond.BOHR, 12AtomBond.BOHR, 15AtomBond.BOHR
        abc = "celldm(1)=10, celldm(2)=1.2, celldm(3)=1.5"
        for (sys, sg, vol) in [("ibrav=5, celldm(1)=10, celldm(4)=0.3", 166, a^3 * sqrt(1 - 3 * 0.3^2 + 2 * 0.3^3)),
                               ("ibrav=6, celldm(1)=10, celldm(3)=1.5", 123, a^2 * c),
                               ("ibrav=7, celldm(1)=10, celldm(3)=1.5", 139, a^2 * c / 2),
                               ("ibrav=8, $abc", 47, a * b * c),
                               ("ibrav=9, $abc", 65, a * b * c / 2),
                               ("ibrav=10, $abc", 69, a * b * c / 4),
                               ("ibrav=11, $abc", 71, a * b * c / 2),
                               ("ibrav=12, $abc, celldm(4)=0.2", 10, a * b * c * sqrt(1 - 0.2^2)),
                               ("ibrav=-12, $abc, celldm(5)=0.2", 10, a * b * c * sqrt(1 - 0.2^2))]
            x = qe(sys)
            @test (sys, spacegroup(x)[1]) == (sys, sg)
            @test cell_volume(x.lattice) ≈ vol
        end
        @test lattice_parameters(qe("ibrav=12, $abc, celldm(4)=0.2").lattice)[6] ≈ acosd(0.2)
        @test lattice_parameters(qe("ibrav=-12, $abc, celldm(5)=0.2").lattice)[5] ≈ acosd(0.2)
        @test lattice_parameters(qe("ibrav=12, A=4, B=5, C=6, cosAB=0.2").lattice)[6] ≈ acosd(0.2)

        # pw.x output: cell in alat or Å, positions in Å or Bohr
        pwo(cell, pos) = read_qe_output(IOBuffer("     Program PWSCF v.7.2\n" *
            "     lattice parameter (alat)  =      10.0000  a.u.\nBegin final coordinates\n$cell\n$pos\nEnd final coordinates\n"))
        q = pwo("CELL_PARAMETERS (alat= 10.00000000)\n 1 0 0\n 0 1 0\n 0 0 1",
                "ATOMIC_POSITIONS (angstrom)\nFe 0 0 0\nFe $(a / 2) $(a / 2) $(a / 2)")
        @test cell_volume(q.lattice) ≈ a^3 && q.frac[2] ≈ [0.5, 0.5, 0.5]
        q = pwo("CELL_PARAMETERS (angstrom)\n 4 0 0\n 0 4 0\n 0 0 4", "ATOMIC_POSITIONS (bohr)\nFe $(2 / AtomBond.BOHR) 0 0")
        @test cell_volume(q.lattice) ≈ 64 && q.frac[1] ≈ [0.5, 0, 0]

        # CIF symmetry from the space-group number or symbol alone
        cif(sym, atoms="Na1 Na 0 0 0\nCl1 Cl 0.5 0.5 0.5") = read_cif(IOBuffer(
            "header before the first block\ndata_x\n_cell_length_a 5.64\n_cell_length_b 5.64\n_cell_length_c 5.64\n" *
            "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n$sym\nloop_\n_atom_site_label\n" *
            "_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n$atoms\n"))
        x = cif("_space_group_IT_number 225")
        @test formula(x) == "Na4Cl4" && spacegroup(x)[1] == 225
        x = @test_logs (:warn, r"several settings") match_mode=:any cif("_space_group_IT_number 227", "Si1 Si 0 0 0")
        @test natoms(x) > 1
        x = @test_logs (:warn, r"several settings") match_mode=:any cif("_symmetry_space_group_name_H-M 'F d -3 m'", "Si1 Si 0 0 0")
        @test natoms(x) > 1
        @test_throws ErrorException cif("_symmetry_space_group_name_H-M 'F d -3 m:9'", "Si1 Si 0 0 0")
        x = @test_logs (:warn, r"no symmetry") match_mode=:any cif("")
        @test formula(x) == "Na1Cl1"
        @test_logs (:warn, r"share a position") match_mode=:any cif("", "Na1 Na 0 0 0\nNa2 Na 0 0 0")

        # VASP 4 POSCAR: no species line, the names come from the comment line
        p = read_poscar(IOBuffer("Si Ge\n1.0\n5.4 0 0\n0 5.4 0\n0 0 5.4\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5\n"))
        @test p.species == [:Si, :Ge]
        @test_throws ErrorException read_poscar(IOBuffer("x\n1.0\n5.4 0 0\n0 5.4 0\n0 0 5.4\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5\n"))

        # XSF: the conventional cell and data grids are skipped
        x = read_xsf(IOBuffer("CRYSTAL\nPRIMVEC\n0 2.7 2.7\n2.7 0 2.7\n2.7 2.7 0\nCONVVEC\n5.4 0 0\n0 5.4 0\n0 0 5.4\n" *
                              "PRIMCOORD\n2 1\nSi 0 0 0\nSi 1.35 1.35 1.35\nCONVCOORD\n8 1\n" * "Si 0 0 0\n"^8 *
                              "BEGIN_BLOCK_DATAGRID_3D\ndensity\nEND_BLOCK_DATAGRID_3D\n"))
        @test natoms(x) == 2 && spacegroup(x)[1] == 227

        # SIESTA cell from LatticeParameters (diamond in its 60° primitive cell)
        s = read_siesta(IOBuffer("LatticeConstant 3.84 Ang\n%block LatticeParameters\n1 1 1 60 60 60\n%endblock LatticeParameters\n" *
                                 "%block ChemicalSpeciesLabel\n1 14 Si\n%endblock ChemicalSpeciesLabel\n" *
                                 "AtomicCoordinatesFormat Fractional\n%block AtomicCoordinatesAndAtomicSpecies\n" *
                                 "0 0 0 1\n0.25 0.25 0.25 1\n%endblock AtomicCoordinatesAndAtomicSpecies\n"))
        @test spacegroup(s)[1] == 227 && cell_volume(s.lattice) ≈ 3.84^3 / sqrt(2)

        # vasprun.xml of an unfinished run: the last structure written
        st(x) = "<structure><crystal><varray name=\"basis\"><v>2.8 0 0</v><v>0 2.8 0</v><v>0 0 2.8</v></varray>" *
                "</crystal><varray name=\"positions\"><v>$x 0 0</v></varray></structure>\n"
        r = read_vasprun(IOBuffer("<modeling><atominfo><array name=\"atoms\"><set><rc><c>Fe</c><c>1</c></rc></set></array>" *
                                  "</atominfo>\n" * st(0) * st(0.1) * "</modeling>"))
        @test r.frac[1] ≈ [0.1, 0, 0] && endswith(r.title, "(last step, unfinished)")

        # .in files: pw.x, FHI-aims or ABINIT by content; unknown files and extensions
        dir = mktempdir()
        f(name, text) = (p = joinpath(dir, name); write(p, text); p)
        @test formula(load_structure(f("aims.in", "lattice_vector 3 0 0\nlattice_vector 0 3 0\nlattice_vector 0 0 3\natom_frac 0 0 0 Fe\n"))) == "Fe1"
        @test formula(load_structure(f("abinit.in", "acell 3*3.0 angstrom\nnatom 1 ntypat 1 typat 1 znucl 26\nxred 0 0 0\n"))) == "Fe1"
        @test_throws ArgumentError load_structure(f("other.in", "nothing to read\n"))
        @test_throws ArgumentError save_structure(joinpath(dir, "x.unknown"), x)
    end

    @testset "printing in the REPL" begin
        si = load_structure(joinpath(STRUCTS, "Si_conv.vasp"))
        @test sprint(show, si) == "Crystal(Si8, 8 atoms)"
        txt = sprint(show, MIME"text/plain"(), supercell(si, 2, 1, 1))
        @test startswith(txt, "Crystal Si16") && occursin("Lattice(a=10.86", txt) && endswith(txt, "… 6 more atoms")
        @test sprint(show, element(:Si)) == "Element(Si, Z=14)"
        @test sprint(show, library_entry("Silicon")) == "LibraryEntry(\"Silicon\", crystal)"
    end

    @testset "view options" begin
        si = load_structure(joinpath(STRUCTS, "Si_conv.vasp"))
        png = joinpath(mktempdir(), "si.png")
        @test snapshot(png, si; settings=ViewSettings(show_labels=true, show_axis_labels=true)) == png && filesize(png) > 0
        v = visualize(si; show=false, panel=false, settings=ViewSettings(style=:line, orthographic=true))
        @test any(p -> p isa AtomBond.Makie.LineSegments, v.plots) && isempty(AtomBond.drawn_atoms(v.display, v.settings))
        cam = AtomBond.camera(v.lscene)
        for (which, dir) in [(:c, [0, 0, 1]), ([1, 1, 0], [1, 1, 0])]
            view_along!(v, which)
            @test norm(cross(normalize(cam.eyeposition[] - cam.lookat[]), normalize(dir))) < 1e-5
        end
        # Turning the view in the screen plane keeps the line of sight; +90° is counter-clockwise.
        view_along!(v, :c)
        d0 = normalize(cam.eyeposition[] - cam.lookat[])
        screen(p) = (d = normalize(cam.eyeposition[] - cam.lookat[]); up = normalize(cam.upvector[]);
                     q = p - cam.lookat[]; (dot(q, cross(up, d)), dot(q, up)))          # (right, up) on screen
        p = cam.lookat[] + 2 * cross(normalize(cam.upvector[]), d0)                    # 2 Å right of centre
        x0, y0 = screen(p)
        rotate_view!(v, 90)
        @test norm(cross(normalize(cam.eyeposition[] - cam.lookat[]), d0)) < 1e-5     # same direction
        x1, y1 = screen(p)
        @test isapprox(x1, 0; atol=1e-4) && isapprox(y1, x0; atol=1e-4) && x0 > 1   # moved from the right to the top
        rotate_view!(v, -90)
        @test all(isapprox.(screen(p), (x0, y0); atol=1e-4))
        rotate_view!(v, 30); rotate_view!(v, 60)
        @test all(isapprox.(screen(p), (x1, y1); atol=1e-4))                         # angles add up
        @test "Turn view 90° clockwise" in [it.label for it in AtomBond.menubar_items(v, "View")]
        # Fewer sphere/bond segments for large structures while viewing; images use full detail
        vs = ViewSettings()
        @test AtomBond.mesh_detail(vs, 100) == (48, 24) && AtomBond.mesh_detail(vs, 1000) == (24, 12) &&
              AtomBond.mesh_detail(vs, 5000) == (16, 8)
        @test AtomBond.mesh_detail(ViewSettings(adaptive_detail=false), 5000) == (48, 24)
        @test AtomBond.mesh_detail(ViewSettings(sphere_resolution=12), 1000) == (12, 12)          # never more than asked
        big = visualize(supercell(si, 4, 4, 4); show=false, panel=false)           # 512 atoms (+ boundary copies)
        @test AtomBond.mesh_detail(big.settings, length(big.display.atoms)) == (24, 12)
        nplots = length(big.plots)
        img = AtomBond.render_image(big)                                         # drawn at full detail, then back
        @test size(img, 1) > 0 && big.settings.adaptive_detail && length(big.plots) == nplots
    end
end
