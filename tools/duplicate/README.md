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
| `--zone-mat <axis> <lo> <hi> <cs> <mat>` | repeatable: links whose bar node has its `<axis>` (x/y/z) coordinate in `[lo,hi]` get `<cs> <mat>` instead; first match wins. Use it for a debonded portion (e.g. the loaded-end sleeve of a lap splice): the link stays, so the bar remains laterally constrained, but it points at a zero-bond material. | none |

Without the new options (`--2d`, `--diameter`, `--bond`) a 3D input gives the
same output as before.

Each link gets a tributary bond length: the mean length of the adjacent bar
elements, halved at a bar end. Its direction vector is the mean direction of
the adjacent bar elements.

## 3D mode

- Bar elements: `truss3d`, `latticeframe3d`. Continuum: `ltrspace`.
- Duplicate: `hangingnode ... dofType 6 2 2 2 2 2 2 dofidmask 6 1 2 3 4 5 6`.
- Link: `latticelink3d` (material e.g. `latticeslip`).
- Writes `Set 1` (bars), `Set 2` (links), `Set 3` (other elements) and the
  fixed sets `Set 4 nodes 1 9`, `Set 5 nodes 1 12`.

## 2D mode

For plane stress / plane strain inputs as written by `t3d2oofem`: header
record with `ndofman`/`nelem`, then the node records, then the element
records, then everything else. That order is kept; no sets are written.

- Bar elements: `truss2d` on ordinary `node` records. Every other record is
  copied unchanged and in place.
- Bar node: rewritten as `hangingnode ... doftype 2 0 2`, i.e. `u` free (slip
  along the bar) and `v` slaved to the host continuum element (the bar follows
  the concrete transversally). A bar node that is already a `hangingnode` is
  left as it is.
- Duplicate: `hangingnode ... doftype 2 2 2` with the bar node's coordinates,
  inserted after the last node record.
- Link: `bondlink1d` (material e.g. `linkslip`, cross section `interfacecs`),
  inserted after the last element record. `ndofman` and `nelem` are updated.

`bondlink1d` couples the `u` DOF only, so all bars must be parallel to x; the
tool stops with an error otherwise. Inclined bars or stirrups would need a 2D
link element.

Example (oofem-examples, `cdpm2-bondslip-2d`):

    duplicate oofemTemp.in oofem.in --diameter 0.052 --bond 3 4

## Compile

    g++ -std=c++17 -O2 -o duplicate duplicate.C
    # plain g++ may be missing on macOS; clang++ works:
    /usr/bin/clang++ -std=c++17 -O2 -o duplicate duplicate.C
