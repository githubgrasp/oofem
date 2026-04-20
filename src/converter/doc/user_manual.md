# Converter User Manual

## Running the converter tests

Tests are diff-based: the converter produces a deterministic OOFEM `.in` file which is compared against a committed `reference.in`.

Build the converter first (from your build directory):

```bash
cd ~/build/oofem
make converter.exec
```

Run all converter tests:

```bash
ctest -R test_converter
```

Run a single test by name:

```bash
ctest -R test_converter_3DCylinderQhull -V
```

The `-V` flag prints the full output, which is useful for diagnosing failures.

### Test layout

Test directories are named `<Case><Mesher>` so the mesher path is obvious from the directory name (e.g. `2DPlateT3D`, `3DFPZQhull`). Each test lives in `src/converter/tests/<TestName>/` and contains:

| File | Purpose |
|------|---------|
| `control.in` | Control file template (also the OOFEM input template) |
| `reference.in` | Committed golden output to diff against |
| `mesh.nodes` + `mesh.voronoi` | Qhull mesh input (qhull path) |
| `mesh.t3d` | T3D mesh input (t3d path) |
| `mesh.sh` | Shell script showing how the mesh was generated (informational) |

### Updating a reference file

After intentionally changing converter output, regenerate the reference from the test directory:

```bash
cd src/converter/tests/3DCylinderQhull
~/build/oofem/src/converter/converter.exec control.in mesh.nodes mesh.voronoi
cp oofem.in reference.in
```

Then re-run the test to confirm it passes.

## Running the converter

The converter needs to know which mesher produced the input. Specify it with
either the `--mesher` (or `-m`) CLI flag, or with a `#@mesher t3d|qhull`
directive in `control.in`. The CLI flag wins if both are present.

```bash
# T3D mesh
converter.exec --mesher t3d control.in mesh.t3d

# Qhull mesh (nodes + voronoi)
converter.exec --mesher qhull control.in mesh.nodes mesh.voronoi

# Qhull mesh (nodes + delaunay + voronoi) — legacy 3-file form
converter.exec --mesher qhull control.in mesh.nodes mesh.delaunay mesh.voronoi
```

If `control.in` already declares the mesher (e.g. all the converter tests do),
the flag can be omitted:

```bash
converter.exec control.in mesh.t3d
converter.exec control.in mesh.nodes mesh.voronoi
```

Output is always written to `oofem.in` in the working directory.

## Control file `#@` directives

The control file is an OOFEM `.in` template. The converter reads it line by line. Lines that start with `#@` are directives — they either configure the converter or mark an injection point for mesh-derived content, and are **not** written to the output. All other lines (including `#` comments and `#%` regression-check blocks) are passed through verbatim.

The available directives differ by path. Use the directives in the table matching the mesh format you are feeding the converter.

### Shared directives (both paths)

These are injection markers, replaced in the output by content derived from the mesh.

| Directive | Purpose |
|-----------|---------|
| `#@INSERT_SETS` | Replaced with `set` records for node/element groups referenced by BCs, loads, and export modules |
| `#@INSERT_LIVELOADS` | Replaced with `NodalLoad`/equivalent records generated from `#@LOAD` requests |
| `#@INSERT_CROSSSECTION` | Replaced with cross-section records for the active grid type |

The counts line (`ndofman … nelem … nset …`) is patched in place once the mesh is processed; it is not an `#@` directive but behaves like one.

### Qhull path directives

Used when the converter is invoked with `mesh.nodes` (+ optional `mesh.delaunay`) and `mesh.voronoi`. Geometry directives set up the region and localiser; they are consumed, not written out.

| Directive | Arguments | Purpose |
|-----------|-----------|---------|
| `#@grid` | `<type>` | Selects the output generator. The unified types are `3dSM` (covers plain, periodic, and fibre-bearing SM) and `3dTM` (mass transport, plain or periodic). A staggered SMTM coupled analysis is expressed as three separate templates — one per subproblem (SM, TM, and the `StaggeredProblem` control file) — each with its own `#@grid` directive. See `Grid::resolveGridType`. |
| `#@mesher` | `t3d` \| `qhull` | Declares which mesher produced the input. Lets `main.C` dispatch without relying on the CLI `--mesher` flag. |
| `#@diam` | `<d>` | Nominal grain diameter — must match the `diam` used when the mesh was generated. |
| `#@perflag` | `3 <px> <py> <pz>` | Periodicity flags per axis (0 = non-periodic, 1 = periodic). For `3dSM`, any axis = 1 switches the writer into periodic mode (control node, `lattice3Dboundary` for boundary-crossing elements). |
| `#@ranint` | `<seed>` | Random integer seed. Non-negative values are replaced with `-time(NULL)`. |
| `#@prism` | `<id> box 6 <xmin> <ymin> <zmin> <xmax> <ymax> <zmax> [refine <r>] [edgerefine <re>] [surfacerefine <rs>] [regionrefine <rr>]` | Defines a box-shaped region. The `*refine` tokens exist for mirroring the generator's `mesh.in` — the converter ignores them, but keeping them in sync keeps the two files readable as a pair. |
| `#@cylinder` | `<id> line 6 <x1> <y1> <z1> <x2> <y2> <z2> …` | Defines a cylindrical region (see `Grid::readQhullControlRecords`). |
| `#@interfacecylinder` | `<id> line 6 <x1> <y1> <z1> <x2> <y2> <z2> …` | Defines a cylindrical interface region. |
| `#@fibre` | `<id> endpoints 6 <x1> <y1> <z1> <x2> <y2> <z2> diameter <d>` | Declares a straight fibre. The converter discretises it into reinforcement nodes at intersections with matrix Voronoi cells, builds `lattice3D` segments along the fibre, and adds `latticelink3D` couplings to the surrounding matrix vertices. |
| `#@notch` | `<id> box 6 <xmin> <ymin> <zmin> <xmax> <ymax> <zmax> material <m>` | Axis-aligned box describing a notch slit. Matrix `lattice3D` / `lattice3Dboundary` elements whose midpoint falls inside the box are written with `crossSect <m> mat <m>` instead of the default `1 1`. Multiple `#@notch` directives may be given; the first matching one wins. The control file is expected to declare a matching `latticecs <m>` / material `<m>` pair. |
| `#@sphereinclusion` | `<id> centre 3 <x> <y> <z> radius <r> itz <t> inside <mi> interface <mif>` | Spherical inclusion with ITZ halo. Effective radius is `r + t/2`. An element whose endpoints both fall inside gets `crossSect mi mat mi`; an element whose endpoints straddle the sphere boundary (one inside, one outside — the ITZ zone) gets `crossSect mif mat mif`. Applied after `#@notch` on the inherited default material. The control file is expected to declare matching `latticecs` / material pairs for both `mi` and `mif`. Use `#@bodyload` to attach a body load to any of the resulting materials. |
| `#@cylinderinclusion` | `<id> line 6 <x1> <y1> <z1> <x2> <y2> <z2> radius <r> itz <t> inside <mi> interface <mif>` | Straight-axis cylindrical inclusion with ITZ halo. The classification uses perpendicular distance from each endpoint to the infinite axis through the two line points; inclusion semantics match `#@sphereinclusion`. Use-case: rebar + ITZ + matrix partitioning for corrosion-cracking studies. |
| `#@bodyload` | `<mat> <bc_id>` | Per-material body load. Any element whose `crossSect` / `mat` resolves to `<mat>` gets `bodyloads 1 <bc_id>` appended. Decoupled from inclusion directives so each test can place its `StructTemperatureLoad` (or similar BC) on the material it needs — e.g. on the matrix for an eigenstrain load (Wong), or on the interface for an eigendisplacement load (corrosion cylinder). Multiple entries allowed (one per material). |
| `#@couplingflag` | *(no arguments)* | Toggle. When present, each emitted element gets `couplingflag 1 couplingnumber <N> <ids…>` appended after the polycoords block. For SM lattice3D / lattice3Dboundary (at Delaunay lines), the ids are the Voronoi-line cross-section elements; for TM latticemt3D (at Voronoi lines), they are the Delaunay-line cross-section elements. Image-side entries (`outsideFlag == 1`) are swapped for their periodic partner via `givePeriodicElement()`. Used by staggered SMTM analyses (e.g. Wong percolation) where SM and TM subproblems exchange per-element data. |
| `#@controlvertex` | `<id> coords 3 <x> <y> <z>` | Declares a specific mesh-node location whose nearest Delaunay-vertex id is exposed at write time via the inline placeholder `#@CTL<id>`. Must also appear as a `controlvertex` entry in the generator's `mesh.in` so the mesher seeds a node at that coordinate. Typical use: naming support/load/monitor nodes whose ids feed into `hpc`, `OutputManager`, BC sets, and check rules. |
| `#@CTL<id>` | *(inline placeholder)* | Substituted on each non-directive line with the mapped OOFEM node id (compact in non-periodic mode, raw in periodic mode). Longer ids substitute first so that e.g. `#@CTL12` isn't shadowed by `#@CTL1`. |
| `#@CTLNODE` | *(inline placeholder, not a line directive)* | Substituted with the numeric id of the periodic control node before each non-directive line is written. Use anywhere the OOFEM template needs that id (e.g. `hpc 2 #@CTLNODE 2 hpcw 1 1.`, `OutputManager tstep_all dofman_output {#@CTLNODE}`, `#NODE number #@CTLNODE dof 2 unknown d`). Inert when `perflag` is fully non-periodic. |
| `#@pov` | *(no arguments)* | Opt in to writing the auxiliary POV-Ray rendering files (`*.vor.line.pov`, `*.vor.cross.pov`, `*.del.line.pov`, `*.del.cross.pov`). Default off — POV files are not generated unless this directive is present. |
| `#@vtk` | *(no arguments)* | Opt in to writing the ParaView `.vtu` files (`*.voronoielement.vtu`, `*.delaunayelement.*.vtu`, `*.fibre.beamElement.vtu`, etc.). Default off — VTK files are not generated unless this directive is present. |

### T3D path directives

Used when the converter is invoked with `mesh.t3d`. These directives describe BCs, loads, and section data attached to named T3D entities (vertices, curves, surfaces).

| Directive | Arguments | Purpose |
|-----------|-----------|---------|
| `#@BC` | `<entType> <entID> <dofMask> <values…>` | Dirichlet BC on a T3D entity. Generates `BoundaryCondition` records and a matching set. |
| `#@LOAD` | `<entType> <entID> <q> [ltf <ltfID>]` | Distributed load on a T3D entity. Emitted via `#@INSERT_LIVELOADS`. |
| `#@DIR` | `<entID> <dx> <dy> <dz>` | Direction vector associated with an entity (used for shell/beam orientation). |
| `#@THICKNESS` | `<entID> <t>` | Shell/plate thickness for an entity. |
| `#@3DSECTION` | `<args>` | 3D section data for shell/beam elements. |
| `#@SHELLWIDTHSCALE` | `<args>` | Per-edge width scaling for shell elements. |
| `#@element` | `<entityKind> <entityID> <elementName> <crossSect> <mat>` | Per-region element override for the T3D writer. `<entityKind>` is one of `vertex`, `curve`, `surface`, `patch`, `shell`. Edges classified to the given entity emit `<elementName>` instead of the writer's hardcoded default (`lattice3D`, `lattice3d`, or `lattice3Dnl` depending on geometry). Edges that don't match any directive use the default. Resolution priority for an edge that touches multiple entities: curve > surface > region. |

`<entType>` uses the numeric codes from `Grid::entityTypeFromString`: 1 vertex, 2 curve, 3 surface, 5 patch, 6 shell.

### Unified `3dSM` writer — conventions for matrix, fibres, and links

`#@grid 3dSM` dispatches to a single writer (`Grid::give3DSMOutput`) that handles plain SM, periodic SM, and periodic SM with fibres, driven entirely by `#@perflag` and the presence of `#@fibre` records. The previously separate `3dFPZ`, `3dFPZFibre`, `3dFibreBenchmark`, `3dGopSha`, `3dTension`, `3dSphere`, and `3dCylinder` grid types have been retired; their analyses are expressed entirely as `control.in` templates using `#@notch` / `#@sphereinclusion` / `#@cylinderinclusion` for per-region material overrides and `#@controlvertex` / `#@CTL<id>` for support/load points.

For the fibre case the writer emits three classes of element, all of them `lattice3D` family, distinguished by their cross-section reference. The control file is expected to declare three `latticecs` records and three materials in this exact slot order:

| Slot | Cross-section / material | Element kind | Geometry source |
|------|--------------------------|--------------|-----------------|
| 1    | `latticecs 1 material 1 shape 2`<br>`latticedamage 1 …` (or any matrix material) | matrix `lattice3D` / `lattice3Dboundary` | polygonal `polycoords` from qhull Voronoi facet |
| 2    | `latticecs 2 material 2 shape 1 radius <r>`<br>`latticelinearelastic 2 …` (or any fibre material) | fibre segment `lattice3D` / `lattice3Dboundary` (Timoshenko frame) | circular cross-section, `radius` taken from `latticecs` |
| 3    | `latticecs 3 material 3`<br>`latticeslip 3 …` (or any bond material) | `latticelink3D` / `latticelink3Dboundary` (fibre↔matrix coupling) | `length`, `diameter`, `dirvector`, `L_end` written by the converter |

The counts header in `control.in` must read `ncrosssect 3 nmat 3 …` whenever fibres are present (drop to 1 for plain SM). The writer prepends `ndofman` and `nelem` automatically.

Reinforcement-node ids are offset by the raw Delaunay-vertex count `N_DelV`: the i-th reinforcement node is written as `node (N_DelV + i)`. The periodic control node id is `N_DelV + N_reinf + 1`, available throughout the template via the `#@CTLNODE` placeholder.

### Unified `3dTM` writer — plain and periodic

`#@grid 3dTM` dispatches to `Grid::give3DTMOutput`, which handles plain and periodic mass transport on the Voronoi dual mesh. When `#@perflag` has any periodic axis, the writer pins the first emitted Voronoi vertex (`bc 1 1`), emits `latticemt3Dboundary` for lines whose endpoints straddle a periodic face (with image-side endpoints swapped for their mirror partners and a `location 2 …` suffix), and appends a control node with `ndofs 3 dofIDmask 3 1 2 3 bc 3 1 1 2` — available in the template via `#@CTLNODE`. Non-periodic behaviour is unchanged (compact 1..N vertex ids, no control node, no boundary elements). Per-element material is resolved through `#@notch` / `#@sphereinclusion` / `#@cylinderinclusion` just like the SM writer, and `#@bodyload` / `#@couplingflag` apply symmetrically.

### Adding a new analysis type

The long-term goal is that supporting a new analysis type means writing a new `control.in`, not touching C++. When adding a directive, register it in the matching reader (`readControlRecords` for T3D, `readQhullControlRecords` for qhull) and add it to the table above. Keep `mesh.in` and the `#@prism`/`#@cylinder`/`#@fibre` etc. lines in `control.in` mirrored so that readers can see the mesh parameters in both files. Prefer extending `give3DSMOutput` / `give3DTMOutput` over adding a new grid type.
