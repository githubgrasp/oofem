Insert node-to-node contact (latticecontact3d) into an existing OOFEM continuum
mesh, either between two meshed bodies or at a free support surface.

Two modes:

Internal interface (between meshed bodies, told apart by cross-section id =
OOFEM region number). An interface = tet faces shared by two tets of different
bodies. The tool duplicates the shared interface nodes (real nodes), makes one
copy per incident body so a node where k bodies meet is split k ways, reassigns
each body's tets to its own copy (unzip; the body listed first keeps the original
id), and inserts one latticecontact3d per adjacent body-pair per node with the
pair's tributary area and area-weighted normal. Every internal face separates
exactly one pair, so tributary area at a triple line is partitioned with no
double count.
  interface in.in out.in --bodies c1,c2,... --contact-cs C --contact-mat M \
            [--pair P:Q[:cs:mat] ...] [--no-set]

  With >=3 bodies a node on a triple line (e.g. two tied lap-spliced bars meeting
  the concrete) splits three ways. Without --pair, every inter-body pair gets a
  contact with the default --contact-cs/--contact-mat. Passing one or more --pair
  restricts contacts to the listed pairs (nodes are still split by full
  incidence), and P:Q:cs:mat overrides the cross-section/material for that pair.
  Example (bar-concrete bond only, no direct bar-bar link):
    interface in.in out.in --bodies 1,2,3 --contact-cs 9 --contact-mat 9 \
              --pair 1:2 --pair 1:3

  Legacy 2-body form (equivalent to --bodies A,B):
  interface in.in out.in --sideA-cs A --sideB-cs B --contact-cs C --contact-mat M [--no-set]

Free-surface interface (specimen against an unmeshed support, e.g. a column
base). The support surface = boundary faces (belonging to one tet) restricted to
a plane or a node set. The tool duplicates those nodes into fixed support nodes,
inserts contacts with the outward normal, and by default TRANSFERS the specimen
nodes' Dirichlet BC to the support duplicates (freeing the specimen side). No tet
unzip. Emits a contact-element set and a support-node set.
  interface in.in out.in --free (--plane <x|y|z> <coord> [--plane-tol t] | --node-set <id>) \
            --contact-cs C --contact-mat M [--no-transfer-bc] [--no-set]
  --no-transfer-bc keeps the specimen BC in place and leaves the support
  duplicates unconstrained (constrain them yourself via the support-node set:
  fixed, prescribed settlement, elastic foundation, ...).

Preconditions:
  * Internal: conforming mesh across the interfaces (adjacent bodies share nodes).
  * The contact LatticeCS and material (e.g. latticebondplast) already exist in
    the mesh; pass their ids. The tool adds only nodes, contacts and sets.

Compile:
  g++ -std=c++17 -O2 -o interface interface.C
  # plain g++ may be missing on macOS; clang++ works:
  /usr/bin/clang++ -std=c++17 -O2 -o interface interface.C
