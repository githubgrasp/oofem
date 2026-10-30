#!/usr/bin/env python3
# Wire loading/support + rigid master-slave sections from LABELS placed in
# oofem.t3d.ctrl. t3d2oofem emits each label as a trailing token on the node line;
# this script resolves label -> node number, so NOTHING is keyed on coordinates or
# dimensions -- change the mesh size or the specimen and it still works. Runs AFTER
# interface (bars = cross-section 2).
#
# Convention (edit the nodeprop labels in oofem.t3d.ctrl to change what's driven):
#   master vertex : m<Role><Tag>   e.g. mSupA, mLoadB     (Role = Sup | Load)
#   slave  cap    : s<Role><Tag>   e.g. sSupA, sLoadB     (same <Role><Tag>)
# Sup masters are fixed (bc 3 1 1 1); Load masters are pulled in z (bc 3 1 1 2).
# Indirect control (analysis record has `hpcmode`): Load masters get bc 3 1 1 0
# load 1 2 (NodalLoad 2 in z) and hpc/hpcw are added to the analysis record so
# that the loaded-end slip of slipLoadA, u_z(bar) - u_z(concrete) at the concrete
# face (the LVDT), is controlled (both are real nodes here, the two copies made
# by the interface tool). The unloaded-end slip is not used: with the stiff bond
# contacts and stiffMode 2 its elastic predictor is ~0 and CALM diverges.
# A metastep record with the token `endcontrol` (removed here) controls the mean
# u_z of the bar nodes instead (more robust at the start, where the slip is ~0).
# Each slave is tied to its cap master via mastermask/doftype (rigid plate, so the
# master carries the whole section reaction). Labels are then stripped. The CHECK
# block is written with the resolved ids: DOFMAN = a Load master (disp), REACTION =
# each Sup master (force).
import re, sys, shutil

f = sys.argv[1] if len(sys.argv) > 1 else "oofem.in"
shutil.copyfile(f, "oofem.nomaster.in")
lines = open(f).read().splitlines()
indirect = any("hpcmode" in ln.split() for ln in lines[:4])

LAB  = re.compile(r'^([ms])(Sup|Load)([A-Za-z0-9]+)$')
SLIP = re.compile(r'^slip([A-Za-z0-9]+)$')            # slip readout: bar vs concrete

# cross-section membership (steel = cs2, concrete = cs1), to pick the copy of a split node
body = {}
for ln in lines:
    if ln.lower().startswith("ltrspace"):
        p = ln.split(); ni = p.index("nodes"); nn = int(p[ni+1]); ns = p[ni+2:ni+2+nn]
        cs = int(p[p.index("crossSect")+1])
        for n in ns: body.setdefault(n, set()).add(cs)

# collect labelled nodes: master[tag], slaves[tag] (steel); slip_bar/slip_con[tag]
master, role, slaves = {}, {}, {}
slip_bar, slip_con = {}, {}
for ln in lines:
    if not ln.startswith("node "): continue
    t = ln.split()
    if len(t) < 8: continue
    nid = t[1]; cs = body.get(nid, ())
    for lab in t[7:]:                                # a node may carry several labels
        m = LAB.match(lab)
        if m:
            if 2 not in cs: continue                 # steel side only
            kind, r, tag = m.group(1), m.group(2), m.group(3)
            if kind == 'm': master[(r, tag)] = nid; role[nid] = r
            else:           slaves.setdefault((r, tag), []).append(nid)
            continue
        s = SLIP.match(lab)
        if s:
            tag = s.group(1)
            if   2 in cs: slip_bar[tag] = nid        # bar (steel) copy
            elif 1 in cs: slip_con[tag] = nid        # concrete copy

for k in master:
    if k not in slaves: sys.exit(f"master {k} has no matching slave cap s{k[0]}{k[1]}")

bc_of, mm_of = {}, {}
for k, M in master.items():
    if k[0] == 'Load': bc_of[M] = "1 1 0 load 1 2" if indirect else "1 1 2"
    else:              bc_of[M] = "1 1 1"
    for s in slaves[k]:
        if s != M: mm_of[s] = M
loadMs = [M for k, M in master.items() if k[0] == 'Load']
supMs  = [M for k, M in master.items() if k[0] == 'Sup']

# slip tags in column order: A, LoadA first (original ld.dat layout), others after
tags = sorted(slip_bar, key=lambda t: (t not in ("A", "LoadA"), t))
slipNodes = [n for tag in tags for n in (slip_bar[tag], slip_con.get(tag)) if n]
if indirect and not ("A" in slip_bar and "B" in slip_bar and "LoadA" in slip_bar and "LoadA" in slip_con):
    sys.exit("hpcmode needs slipA, slipB and slipLoadA labels in oofem.t3d.ctrl")

# rewrite node lines: strip the trailing label; masters get bc, slaves mastermask
out = []
for ln in lines:
    if indirect and "hpcmode" in ln.split():
        ln = re.sub(r'\s*hpcw?\s+\d+(\s+[-\d.eE+]+)*', '', ln)   # drop old hpc/hpcw
        if "endcontrol" in ln.split():          # mean u_z of the bar nodes
            ln = re.sub(r'\s*\bendcontrol\b', '', ln)
            out.append(f"{ln} hpc 4 {slip_bar['A']} 3 {slip_bar['B']} 3 hpcw 2 0.5 0.5"); continue
        out.append(f"{ln} hpc 4 {slip_bar['LoadA']} 3 {slip_con['LoadA']} 3 hpcw 2 1 -1"); continue
    if ln.startswith("node "):
        t = ln.split(); nid = t[1]
        base = " ".join(t[:7])                        # node id coords 3 x y z
        if any(LAB.match(x) or SLIP.match(x) for x in t[7:]):
            if   nid in bc_of: out.append(f"{base} bc 3 {bc_of[nid]}"); continue
            elif nid in mm_of: M = mm_of[nid]; out.append(f"{base} mastermask 3 {M} {M} {M} doftype 3 1 1 1"); continue
            else:              out.append(base); continue   # slave concrete copy / slip node -> strip label
        out.append(ln); continue
    if ln.startswith("OutputManager"):
        ln = re.sub(r'\s*dofman_output\s*\{[^}]*\}', '', ln).rstrip()
        ln = f"{ln} dofman_output {{ {' '.join(dict.fromkeys(loadMs + supMs + slipNodes))} }}"
    out.append(ln)

txt = re.sub(r'\n?#%BEGIN_CHECK%.*?#%END_CHECK%\n?', '\n', "\n".join(out), flags=re.S).rstrip() + "\n"
chk = ["#%BEGIN_CHECK%", f"#DOFMAN number {loadMs[0]} dof 3 type d"]
chk += [f"#REACTION number {m} dof 3" for m in supMs]
for tag in tags:                                      # slip = u_z(bar) - u_z(concrete)
    b, c = slip_bar[tag], slip_con.get(tag)
    if c: chk += [f"# slip {tag}: bar {b} minus concrete {c} (dof 3)",
                  f"#DOFMAN number {b} dof 3 type d",
                  f"#DOFMAN number {c} dof 3 type d"]
chk += ["#%END_CHECK%"]
open(f, "w").write(txt + "\n".join(chk) + "\n")

print("resolved masters:", {f"{k[0]}{k[1]}": v for k, v in master.items()})
print(f"CHECK: DOFMAN {loadMs[0]} ; REACTION {supMs}")
print(f"slip readouts (bar,concrete): {[(t, slip_bar[t], slip_con.get(t)) for t in tags]}")
