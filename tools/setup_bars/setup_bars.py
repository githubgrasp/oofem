#!/usr/bin/env python3
# Finish the frame/truss TensileSplice model after `duplicate`. Shared by ../Truss
# and ../Beam (same file in both). Like ../Solid/0p02/setup_masters.py it is
# LABEL-driven: oofem.t3d.ctrl tags bar-end vertices with nodeprop labels, which
# t3d2oofem writes as trailing tokens on the node line (several per node allowed,
# e.g. `nodeprop mLoadA slipLoadA`). Nothing here is keyed on
# coordinates, so it works at any mesh size.
#
# Labels (edit them in oofem.t3d.ctrl):
#   m<Role><Tag>  bar end, Role = Sup (fixed, Set 4) or Load (pulled in z, Set 5)
#   slip<Tag>     bar node for the slip readout; its duplicate (hanging node in
#                 the concrete, found via the link element) is output as well:
#                 slip = u_z(bar) - u_z(concrete)
#   crk<Tag><M|P> concrete surface node below (M) / above (P) a crack; gauge Tag
#                 opening = u_y(P) - u_y(M)
#
# What duplicate leaves and this script fixes:
#   - records it does not know (e.g. Isole, latticeframeelastic) end up after the
#     sets -> all records are put back into OOFEM order
#   - it drops the ctrl's Set records and writes Set 4 = node 9, Set 5 = node 12
#     -> Set 4/5 are rewritten from the labels (Set 1 bars, 2 links, 3 rest kept)
#   - header counts (ndofman, nelem, ncrosssect, nmat, nbc, nltf, nset) recomputed
# Then the labels are stripped, dofman_output and the extractor CHECK block written.
#
# Indirect control: if the analysis record has `hpcmode`, `hpc`/`hpcw` are added
# so that the mean slip of slipA and slipB is controlled: u_z of the bar node
# minus u_z of the nearest concrete (ltrspace) node. The duplicate itself is a
# slave node and cannot be used in hpc (no equation number). With
# the token `crackcontrol` in the record (removed here), the mean opening of the
# crk gauges is controlled instead; `crackcontrol <beta>` controls
# mean u_z + beta * mean opening. Each metastep record is treated separately.
#
# --truss (only for truss bars with latticelink3d links): truss nodes have no
# rotational stiffness, but latticelink3d needs six DOFs per node. duplicate
# embeds only the translations of a truss duplicate (doftype 6 2 2 2 0 0 0), so
# the rotations of the bar nodes and of their duplicates are fixed here
# (bc 6 0 0 0 1 1 1, BC 1 = zero). Not needed with bondlink3dtruss links
# (duplicate --link bondlink3dtruss: translations only) or frame (lattice3d)
# bars (their duplicates are continuumframenodes).
#
# usage: python3 setup_bars.py [--truss] [oofem.in]
import re, sys, shutil

args = [a for a in sys.argv[1:] if a != "--truss"]
truss = "--truss" in sys.argv[1:]
f = args[0] if args else "oofem.in"
shutil.copyfile(f, "oofem.nobars.in")
lines = open(f).read().splitlines()

DOFMAN = {"node", "hangingnode", "continuumframenode"}
ELEM   = {"ltrspace", "truss3d", "lattice3d", "latticelink3d", "bondlink3dtruss"}
CS     = {"simplecs", "latticecs", "interfacecs"}
BC     = {"boundarycondition", "nodalload"}
LTF    = {"constantfunction", "piecewiselinfunction"}
LAB  = re.compile(r'^m(Sup|Load)([A-Za-z0-9]+)$')
SLIP = re.compile(r'^slip([A-Za-z0-9]+)$')
CRK  = re.compile(r'^crk([A-Za-z0-9]+)([MP])$')

# header = everything up to and including the ndofman/nelem record
h = next(i for i, l in enumerate(lines) if "ndofman" in l and "nelem" in l)
head, body = lines[:h + 1], lines[h + 1:]

rec = {k: [] for k in ("dofman", "elem", "cs", "mat", "bc", "ltf", "set", "rest")}
for l in body:
    t = l.split()
    if not t or t[0].startswith("#"):
        rec["rest"].append(l); continue
    k = t[0].lower()
    if   k in DOFMAN: rec["dofman"].append(l)
    elif k in ELEM:   rec["elem"].append(l)
    elif k in CS:     rec["cs"].append(l)
    elif k in BC:     rec["bc"].append(l)
    elif k in LTF:    rec["ltf"].append(l)
    elif k == "set":  rec["set"].append(l)
    else:             rec["mat"].append(l)      # anything else is a material

# bar node -> its hanging duplicate, from the links
dup = {}
for l in rec["elem"]:
    t = l.split()
    if t[0] in ("latticelink3d", "bondlink3dtruss"):
        dup[t[4]] = t[5]
dupNodes = set(dup.values())

# resolve labels on node records, strip them
sup, load, slip, crk = [], [], {}, {}
dofman = []
for l in rec["dofman"]:
    t = l.split()
    labs = [x for x in t[7:] if LAB.match(x) or SLIP.match(x) or CRK.match(x)] if t[0] == "node" else []
    for x in labs:                             # a node may carry several labels
        m, s, c = LAB.match(x), SLIP.match(x), CRK.match(x)
        if m: (sup if m.group(1) == "Sup" else load).append(t[1])
        elif s: slip[s.group(1)] = t[1]
        else: crk.setdefault(c.group(1), {})[c.group(2)] = t[1]
    if labs:
        l = " ".join(t[:7] + [x for x in t[7:] if x not in labs])
    dofman.append(l)
if not sup or not load:
    sys.exit("need at least one mSup* and one mLoad* label in oofem.t3d.ctrl")

if truss:
    ends = set(sup + load)                     # rotations already fixed by Set 4/5
    for i, l in enumerate(dofman):
        t = l.split()
        if t[0] == "hangingnode" and t[1] in dupNodes:
            dofman[i] = l + " bc 6 0 0 0 1 1 1"
        elif t[0] == "node" and t[1] in dup and t[1] not in ends:
            dofman[i] = " ".join(t[:7]) + " dofidmask 6 1 2 3 4 5 6 bc 6 0 0 0 1 1 1"

# nearest concrete (ltrspace) node to each slip bar node, for the slip control
xyz = {t[1]: [float(x) for x in t[4:7]] for t in (l.split() for l in dofman) if "coords" in t[:3]}
conc = {n for t in (l.split() for l in rec["elem"]) if t[0] == "ltrspace"
        for n in t[t.index("nodes") + 2:t.index("nodes") + 6]}
near = {b: min(conc, key=lambda n: sum((p - q) ** 2 for p, q in zip(xyz[n], xyz[b])))
        for b in slip.values()}

# slip tags in column order: A, LoadA first (original ld.dat layout), others after
tags = sorted(slip, key=lambda t: (t not in ("A", "LoadA"), t))
slipNodes = [n for tag in tags for n in (slip[tag], dup.get(slip[tag])) if n]
bad = [g for g in crk if set(crk[g]) != {"M", "P"}]
if bad:
    sys.exit(f"crack gauges without both M and P node: {bad}")
gauges = sorted(crk)
crkNodes = [crk[g][s] for g in gauges for s in "MP"]

sets = [s for s in rec["set"] if s.split()[1] not in ("4", "5")]
sets += [f"Set 4 nodes {len(sup)} {' '.join(sup)}", f"Set 5 nodes {len(load)} {' '.join(load)}"]
sets.sort(key=lambda s: int(s.split()[1]))

cnt = {"ndofman": len(dofman), "nelem": len(rec["elem"]), "ncrosssect": len(rec["cs"]),
       "nmat": len(rec["mat"]), "nbc": len(rec["bc"]), "nic": 0, "nltf": len(rec["ltf"]),
       "nset": len(sets)}
head[-1] = " ".join(f"{k} {v}" for k, v in cnt.items())

out = []
for l in head:
    if "hpcmode" in l.split():
        l = re.sub(r'\s*hpcw?\s+\d+(\s+[-\d.eE+]+)*', '', l)   # drop old hpc/hpcw
        m = re.search(r'\s*\bcrackcontrol\b(\s+[-\d.eE+]+)?', l)
        dofs, w = [], []
        if not m or m.group(1):                # mean slip of the unloaded bar ends
            if not ("A" in slip and "B" in slip):
                sys.exit("hpcmode needs slipA and slipB labels in oofem.t3d.ctrl")
            dofs = [f"{n} 3" for t in "AB" for n in (slip[t], near[slip[t]])]
            w = [0.5, -0.5, 0.5, -0.5]
        if m:                                  # (+ beta *) mean opening (u_y P - M)
            if not gauges:
                sys.exit("crackcontrol needs crk<Tag>M/P labels in oofem.t3d.ctrl")
            l = l[:m.start()] + l[m.end():]
            b = (float(m.group(1)) if m.group(1) else 1.) / len(gauges)
            dofs += [f"{n} 2" for n in crkNodes]
            w += [-b, b] * len(gauges)
        l += f" hpc {2 * len(dofs)} {' '.join(dofs)} hpcw {len(w)} {' '.join(f'{x:g}' for x in w)}"
    if l.startswith("OutputManager"):
        l = re.sub(r'\s*dofman_output\s*\{[^}]*\}', '', l).rstrip()
        l = f"{l} dofman_output {{ {' '.join(dict.fromkeys(load + sup + slipNodes + [near[slip[t]] for t in tags] + crkNodes))} }}"
    out.append(l)
out += dofman + rec["elem"] + rec["cs"] + rec["mat"] + rec["bc"] + rec["ltf"] + sets

chk = ["#%BEGIN_CHECK%", f"#DOFMAN number {load[0]} dof 3 type d"]
chk += [f"#REACTION number {n} dof 3" for n in sup]
for tag in tags:                               # slip = u_z(bar) - u_z(concrete)
    b, c = slip[tag], dup.get(slip[tag])
    if c: chk += [f"# slip {tag}: bar {b} minus concrete {c} (dof 3)",
                  f"#DOFMAN number {b} dof 3 type d",
                  f"#DOFMAN number {c} dof 3 type d"]
for g in gauges:                               # opening = u_y(P) - u_y(M)
    chk += [f"# crack gauge {g}: {crk[g]['P']} minus {crk[g]['M']} (dof 2)",
            f"#DOFMAN number {crk[g]['M']} dof 2 type d",
            f"#DOFMAN number {crk[g]['P']} dof 2 type d"]
chk += ["#%END_CHECK%"]
open(f, "w").write("\n".join(out) + "\n" + "\n".join(chk) + "\n")

print("header:", head[-1])
print("slip control (bar, nearest concrete node, distance mm): " + ", ".join(
    f"{t} {slip[t]} {near[slip[t]]} {1e3 * sum((p - q) ** 2 for p, q in zip(xyz[near[slip[t]]], xyz[slip[t]])) ** .5:.1f}"
    for t in tags))
print(f"Sup {sup}  Load {load}  slip (bar,concrete): {[(t, slip[t], dup.get(slip[t])) for t in tags]}")
if gauges:
    print(f"crack gauges (M,P): {[(g, crk[g]['M'], crk[g]['P']) for g in gauges]}")
