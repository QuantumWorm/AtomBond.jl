# Working from Julia code

Everything the window does is also available as functions. Structures are
immutable [`Crystal`](@ref) values, and editing functions return new ones.
The examples on this page run when the documentation is built, so their
output is real. They use the example files that come with the package:

```@example main
using AtomBond
example(name) = joinpath(pkgdir(AtomBond), "examples", "structures", name)
si = load_structure(example("Si_conv.vasp"))
```

## Structure and symmetry

```@example main
natoms(si), formula(si), lattice_parameters(si.lattice)
```

```@example main
rutile = load_structure(example("TiO2_rutile.cif"))     # symmetry-expanded CIF
spacegroup(rutile)
```

```@example main
natoms(primitive_cell(si)), natoms(conventional_cell(si))
```

```@example main
natoms(supercell(si, 2, 2, 1)), natoms(supercell(si, [1 -1 0; 1 1 0; 0 0 1]))
```

## Bonds and measurements

```@example main
bonds = find_bonds(si)
length(bonds), round(bonds[1].length; digits=4)
```

```@example main
r = cartesian_positions(si)
b = first(filter(b -> b.i == 1, bonds))             # a bond of atom 1
round(distance(r[1], r[b.j] + si.lattice.matrix * b.shift); digits=4)
```

```@example main
water = load_library("Water")
p = cartesian_positions(water)                        # O, H, H
round(bond_angle(p[2], p[1], p[3]); digits=2)         # H–O–H angle
```

```@example main
round(min_image_distance(si, 1, 8); digits=4)         # shortest periodic distance
```

## Editing

```@example main
sige = substitute(si, [5], :Ge)
formula(sige)
```

```@example main
c = delete_atoms(si, [1])
c = move_atoms(c, [4], [0.1, 0, 0])                   # Cartesian Å
c = add_atom(c, :O, [0.5, 0.5, 0.5])                  # fractional
formula(c)
```

```@example main
gaas = load_library("GaAs")
formula(replace_species(gaas, :As, :P))
```

```@example main
gaas2 = set_position(gaas, 5, (0.3, 0.25, 0.25))
gaas2.frac[5]
```

```@example main
shifted = translate_atoms(si, [0.5, 0, 0]; fractional=true)
length(atoms_outside_cell(shifted)), length(atoms_outside_cell(wrap(shifted)))
```

## Surfaces and superlattices

```@example main
s = slab(si, (1, 1, 1); layers=4, vacuum=15)          # termination chosen automatically
natoms(s), s.pbc
```

```@example main
ge = Crystal(cubic(5.658), fill(:Ge, 8), si.frac)
sl = superlattice(si, ge, (0, 0, 1); layers=(4, 4), strain=true)
formula(sl)
```

## Constraints and DFT inputs

```@example main
relax = with_constraints(fix_atoms(s, 1:4); cell=(0, 0, 0, 0, 0, 0))
count(f -> all(f), relax.fixed), describe_cell_constraints(relax.cell_constraints)
```

```@example main
kpoint_grid(si; kspacing=0.25)
```

```@example main
dir = mktempdir()
write_qe(joinpath(dir, "si.in"), si; ecutwfc=40, pseudo="{El}.UPF")
write_abinit(joinpath(dir, "si.abi"), si; ecutwfc=40, pseudo="{El}.psp8")
write_siesta(joinpath(dir, "si.fdf"), si)
readdir(dir)
```

```@example main
save_structure(joinpath(dir, "SiGe.cif"), sige)
formula(load_structure(joinpath(dir, "SiGe.cif")))
```

```@example main
frames = read_castep_geom_frames(example("Si_relax.geom"))   # every step of a relaxation
length(frames), first(frames).title
```

## The window from code

These open a window, so they are not run when the documentation is built:

```julia
v = visualize(si)                               # window with the control panel
view_along!(v, :c); view_along!(v, [1, 1, 0])   # a, b, c, a*, b*, c* or [uvw]
v.settings.bond_stubs = false; redraw!(v)       # any setting, then redraw
save_view("si_view.png", v)                     # just the 3D view
snapshot("si.png", si; settings=ViewSettings(style=:cpk))   # PNG, no window

select!(v, 5); select!(v, 1; add=true)          # same as click / Shift+click
v.info[]                                        # "Si5 – Si1   distance 2.3517 Å"
edit!(c -> substitute(c, [5], :Ge), v, "make SiGe")   # undoable in the window
undo!(v); redo!(v)
set_crystal!(v, supercell(si, 2, 2, 2))         # show another structure, not undoable
new_tab!(v, load_library("GaN"))                # a second structure in a new tab
```
