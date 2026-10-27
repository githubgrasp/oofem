# duplicate: bond links between bars and continuum

Post-processes an OOFEM input file in which reinforcement bars (truss or frame
elements) share the domain with a continuum mesh, but not its nodes. Each bar
node is duplicated. The duplicate becomes a hanging node constrained to the
continuum element it lies in, and the bar node and its duplicate are joined by
a bond link element. The bond-slip law is then set by the link's material.

Typical pipeline:

    t3d -i mesh.in -o mesh.out
    t3d2oofem oofem.t3d.ctrl mesh.out oofemTemp.in
    duplicate oofemTemp.in oofem.in [options]

## Usage

    duplicate inputfile outputfile [--2d | --3d] [--diameter <d>] [--bond <cs> <mat>]
              [--zone-mat <axis> <lo> <hi> <cs> <mat>]...

| option | meaning | default |
|---|---|---|
| `--2d`, `--3d` | force the mode | from the `domain` record: `2d...` or `planestrain` gives 2D, anything else 3D |
| `--diameter <d>` | bond diameter written to every link (bond area = pi d length). For n bars lumped into one truss, use n times the bar diameter. | 0.016 |
| `--bond <cs> <mat>` | cross section and material of the links | 4 4 |
| `--zone-mat <axis> <lo> <hi> <cs> <mat>` | repeatable: bond length in the band `[lo,hi]` of the `<axis>` (x/y/z) coordinate gets `<cs> <mat>` instead; first match wins. The tributary length of a bar node (one half of each adjacent bar element) is assigned half by half, by the midpoint of each half, so a band edge at a bar node is placed exactly there: that node gets one link per material (same duplicate), each with its part of the length. Use it for a debonded portion (e.g. the loaded-end sleeve of a lap splice): the link stays, so the bar remains laterally constrained, but it points at a zero-bond material. | none |


Without the new options (`--2d`, `--diameter`, `--bond`) a 3D input gives the
same output as before.

Each link gets a tributary bond length: the mean length of the adjacent bar
elements, halved at a bar end. Its direction vector is the mean direction of
the adjacent bar elements.

## 3D mode

- Bar elements: `truss3d`, `lattice3d` (frame; replaces the obsolete
  `latticeframe3d`). Continuum: `ltrspace`. All `lattice3d` records are
  treated as bars, so this mode is not meant for lattice models of the
  concrete itself.
- Duplicate of a `lattice3d` node: `continuumframenode ... dofType 6 2 2 2 2 2 2
  dofidmask 6 1 2 3 4 5 6`. Translations and rotations follow the host
  tetrahedron (rotation = 1/2 curl u).
- Duplicate of a `truss3d` node: `hangingnode ... dofType 6 2 2 2 0 0 0
  dofidmask 6 1 2 3 4 5 6`. Only the translations are embedded. The link needs
  six DOFs per node, but a truss node has no rotational stiffness, so the
  rotations of the bar nodes and their duplicates have to be fixed in the
  input (e.g. `bc 6 0 0 0 1 1 1` with a zero boundary condition 1).
- Link: `latticelink3d` (material e.g. `latticeslip`).
- Writes `Set 1` (bars), `Set 2` (links), `Set 3` (other elements) and the
  fixed sets `Set 4 nodes 1 9`, `Set 5 nodes 1 12`.

## 2D mode

For plane stress / plane strain inputs as written by `t3d2oofem`: header
record with `ndofman`/`nelem`, then the node records, then the element
records, then everything else. That order is kept; no sets are written.

- Bar elements: `truss2d`, or `lattice2d` used as a frame (beam) element.
  Every other record is copied unchanged and in place; the bar nodes stay
  ordinary nodes. All `lattice2d` records are treated as bars, so this mode
  is not meant for lattice models of the concrete itself. `lattice2d`
  without `gpCoords` takes the element midpoint as integration point (zero
  eccentricity), so `t3d2oofem` can write it directly with
  `edgetype lattice2d`.
- Duplicate (`truss2d` bars): `hangingnode ... doftype 2 2 2` with the bar
  node's coordinates, inserted after the last node record.
- Link (`truss2d` bars): `bondlink2dtruss` between bar node (first) and
  duplicate (second), inserted after the last element record. The link
  resolves the relative displacement into slip along the bar and a lateral
  component, so bars may have any direction in the x-y plane. Use it with
  e.g. `linkslip` (bond-slip law for the slip, `kl` for the lateral
  stiffness) and an `interfacecs` cross section. `ndofman` and `nelem` are
  updated.

Example (oofem-examples, `cdpm2-bondslip-2d`):

    duplicate oofemTemp.in oofem.in --diameter 0.052 --bond 3 4

The same call is used with `lattice2d` bars (`cdpm2-reinforced-beam-2d`).

### 2D mode with `lattice2d` bars

`lattice2d` bars (frame elements: u, v and rotation per node) get the same
`hangingnode` and `bondlink2dtruss` as `truss2d` bars. The link uses only the
translations of both nodes, so the rotation of the bar node is not tied to
the concrete; it is carried by the frame elements.

## Compile

    g++ -std=c++17 -O2 -o duplicate duplicate.C
    # plain g++ may be missing on macOS; clang++ works:
    /usr/bin/clang++ -std=c++17 -O2 -o duplicate duplicate.C
