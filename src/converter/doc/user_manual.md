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
