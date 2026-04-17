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
ctest -R test_converter_3DCylinder -V
```

The `-V` flag prints the full output, which is useful for diagnosing failures.

### Test layout

Each test lives in `src/converter/tests/<TestName>/` and contains:

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
cd src/converter/tests/3DCylinder
~/build/oofem/src/converter/converter.exec control.in mesh.nodes mesh.voronoi
cp oofem.in reference.in
```

Then re-run the test to confirm it passes.

## Running the converter

```bash
# T3D mesh
converter.exec control.in mesh.t3d

# Qhull mesh (nodes + voronoi)
converter.exec control.in mesh.nodes mesh.voronoi

# Qhull mesh (nodes + delaunay + voronoi)
converter.exec control.in mesh.nodes mesh.delaunay mesh.voronoi
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
| `#@grid` | `<type>` | Selects the output generator (e.g. `3dSM`, `3dTM`, `3dCoupledSMTM`, `3dCylinder`, …). See `Grid::resolveGridType` for the full list. |
| `#@diam` | `<d>` | Nominal grain diameter — must match the `diam` used when the mesh was generated. |
| `#@perflag` | `3 <px> <py> <pz>` | Periodicity flags per axis (0 = non-periodic, 1 = periodic). |
| `#@ranint` | `<seed>` | Random integer seed. Non-negative values are replaced with `-time(NULL)`. |
| `#@prism` | `<id> box 6 <xmin> <ymin> <zmin> <xmax> <ymax> <zmax> [refine <r>] [edgerefine <re>] [surfacerefine <rs>] [regionrefine <rr>]` | Defines a box-shaped region. The `*refine` tokens exist for mirroring the generator's `mesh.in` — the converter ignores them, but keeping them in sync keeps the two files readable as a pair. |
| `#@cylinder` | `<id> line 6 <x1> <y1> <z1> <x2> <y2> <z2> …` | Defines a cylindrical region (see `Grid::readQhullControlRecords`). |
| `#@interfacecylinder` | `<id> line 6 <x1> <y1> <z1> <x2> <y2> <z2> …` | Defines a cylindrical interface region. |

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

`<entType>` uses the numeric codes from `Grid::entityTypeFromString`: 1 vertex, 2 curve, 3 surface, 5 patch, 6 shell.

### Adding a new analysis type

The long-term goal is that supporting a new analysis type means writing a new `control.in`, not touching C++. When adding a directive, register it in the matching reader (`readControlRecords` for T3D, `readQhullControlRecords` for qhull) and add it to the table above. Keep `mesh.in` and the `#@prism`/`#@cylinder` etc. lines in `control.in` mirrored so that readers can see the mesh parameters in both files.
