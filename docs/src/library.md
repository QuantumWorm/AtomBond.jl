# Structure library and element data

## Structure library

File → *Structure library…* (Ctrl+L, also Build → *From library…*) lists
known crystals and molecules with their formula, symmetry, cell and source
(database entry and paper); *Open* loads one in a new tab. The structure is
Untitled until saved, so Save never writes into the library. In code:
[`structure_library`](@ref), [`load_library`](@ref) (by name or formula),
[`open_library!`](@ref).

**Finding a structure:** *Show* picks crystals, molecules or both; *Group*
narrows the list to one group; *Search* keeps the structures whose name,
formula, group or use contain all the words typed (`ethane`, `WF6`,
`halide`, `etch`, `precursor`).

| Group | Structures | Source |
|-------|------------|--------|
| Semiconductors: group IV | Si, diamond, Si0.5Ge0.5 (ordered model) | COD (Hom et al. 1975 precision lattice constants; Dismukes 1964 for SiGe) |
| Semiconductors: III–V | GaAs, AlAs, InAs, GaP, AlP, InP, GaSb, AlSb, InSb (zincblende); GaN, AlN, InN (wurtzite) | 300 K lattice constants from Vurgaftman, Meyer, Ram-Mohan, *J. Appl. Phys.* 89, 5815 (2001); wurtzite u from Paszkowicz et al. (2003, 2004) |
| Oxides | α-quartz SiO₂, corundum Al₂O₃, monoclinic HfO₂ and ZrO₂, rutile and anatase TiO₂, ZnO, MgO, β-Ga₂O₃, SrTiO₃ | COD (published refinements, listed per entry) |
| Carbon: graphite and graphene | graphite (2H); graphene (one layer, periodic in-plane) | COD (Trucano & Chen 1975, neutron diffraction); graphene from that lattice constant |
| Carbon: nanotubes | single-wall (5,5), (10,10), (9,0), (10,0), (6,5) | Rolled graphene sheet with graphite's C–C distance: models, periodic along the tube axis only |
| Molecules | 117 molecules in the groups below | Experimental gas-phase geometries from spectroscopy and electron diffraction, cited per molecule |

Molecules are grouped the way chemistry textbooks classify them: first
organic (carbon compounds), organometallic (metal–carbon bonds) or inorganic,
then by functional group or type of compound. Within a group they are in
alphabetical order.

| Class | Groups (examples) |
|-------|-------------------|
| Organic | hydrocarbons (CH₄, C₂H₂, propyne, cyclopropane, C₆H₆); halocarbons (CF₄, CH₂Cl₂, CHCl₃); alcohols and ethers (CH₃OH, oxirane); carbonyl compounds: aldehydes, ketones, acids, esters, amides (H₂CO, HCOOH, ketene, methyl formate, formamide); nitrogen and sulfur compounds (HCN, CH₃CN, DMSO); heterocycles (pyridine, furan, thiophene); radicals (CH, CH₂, HCO, CF₃) |
| Organometallic | metal alkyls (Al(CH₃)₃ and its dimer, Zn(CH₃)₂, In(CH₃)₃); metal carbonyls (Ni(CO)₄, W(CO)₆) |
| Inorganic | elements (H₂, O₂, O₃, Cl₂, P₄, S₈); hydrides (H₂O, NH₃, SiH₄, B₂H₆, HF); oxides and sulfides (CO, CO₂, SO₂, SO₃, NO₂, N₂O₄, SiO, CS₂); oxoacids (HNO₃, H₂SO₄, HONO, HOF); nonmetal halides (BCl₃, SiCl₄, NF₃, SF₆, OF₂); metal halides (TiCl₄, WF₆, Al₂Cl₆); oxyhalides (COCl₂, POCl₃, SOCl₂, SO₂Cl₂); radicals (NH₂, SiH₂, SiH₃, HO₂, SH, ClO) |
| Carbon | fullerenes (C₆₀) |

Molecules used in semiconductor processing also say what for (*Used as*):
etch gas, dopant source, deposition precursor, ALD/CVD co-reactant
(oxidant, nitrogen source), wet etch and cleaning, or solvent. Search finds them by that too.

- **Crystals** are CIFs from the
  [Crystallography Open Database](https://www.crystallography.net), stored as
  downloaded.
- **III–V semiconductors** come from one review at 300 K instead, because
  COD's III–V entries are mostly old tables and some are off.
- **Molecules** are experimental gas-phase geometries: the bond lengths and
  angles of the cited papers and books (each value read in the authors' own
  text, never taken from a data compilation), from which `scripts/molecules.jl`
  computes the coordinates. Each entry gives its point group and references;
  a few carry a note where the structure is approximate. They open in a
  non-periodic box with 5 Å of space; the
  library's *Molecules: in a periodic cubic box* switch makes the box
  periodic for plane-wave codes.
- **Radicals** are open-shell species (surface adsorbates, plasma and CVD
  intermediates); their note gives the spin state to use.
- **Graphene and nanotubes** are periodic only along some directions
  (extended XYZ with `pbc`): graphene in-plane, a nanotube along its axis
  (c), with vacuum elsewhere. Nanotubes are ideal rolled sheets, so relax them
  before use.
- **SiGe:** the real alloy is random; the library holds the ordered model
  (Si and Ge on the two diamond sublattices) at the measured alloy lattice
  constant.

The library is `data/library/library.toml` with `crystals/` and `molecules/`
beside it; `scripts/build_library.jl` downloads the COD files and rebuilds
it, with the molecules from `scripts/molecules.jl`.
Your own structures go into `~/.atombond/library/library.toml` (or
`$ATOMBOND_HOME/library/`), in the same format with the files next to it (`category` sets the group;
an optional `use` is shown and searched).

## Element data

All element data is in one table, `data/elements.csv` (Z 1–103), generated by
`scripts/build_element_table.jl` from the published sources:

| Column | Used for | Source |
|--------|----------|--------|
| `mass` | QE and LAMMPS masses | NIST *Atomic Weights and Isotopic Compositions* (IUPAC/CIAAW standard atomic weights) |
| `covalent_radius` | bonds, ball size | Cordero et al., *Dalton Trans.* 2008, 2832 |
| `vdw_radius_bondi` | CPK size | Bondi, *J. Phys. Chem.* 68, 441 (1964), with Mantina et al. 2009 |
| `vdw_radius_alvarez` | CPK size where Bondi has none | Alvarez, *Dalton Trans.* 42, 8617 (2013) |
| `color_jmol` | atom colours | Jmol element colours (H white, C grey, N blue, O red, …) |

[`element`](@ref)`(:Si)` gives an [`Element`](@ref);
[`element_property`](@ref)`(:Si, :vdw_radius_alvarez)` reads any column.

**Your own values:** put an `elements.csv` with the same layout in
`~/.atombond/` (or `$ATOMBOND_HOME`), with only the columns and rows to
change:

```
symbol,color,vdw_radius
O,FF2020,1.40
Si,1B3BFA,
```

`color` and `vdw_radius` take precedence over the built-in columns, and new
columns become readable with `element_property`.
[`use_vesta_colors!`](@ref) takes the colours (and optionally the radii) from
a VESTA `elements.ini`.
