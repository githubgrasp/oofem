~/build/oofem/src/generator/generator.exec mesh.in
mv nodes.dat mesh.nodes
qvoronoi p Fv < mesh.nodes > mesh.voronoi

# The registered ctest runs only the SM subproblem (control.in → oofem.in,
# diffed against reference.in). The other two templates (control.tm.in and
# control.smtm.in) are kept as examples showing how the full staggered-SMTM
# pipeline assembles. The SM invocation is run LAST so `oofem.in` ends up
# holding the SM output for the ctest diff.
~/build/oofem/src/converter/converter.exec control.tm.in   mesh.nodes mesh.voronoi && mv oofem.in oofem.tm.in
~/build/oofem/src/converter/converter.exec control.smtm.in mesh.nodes mesh.voronoi && mv oofem.in oofem.smtm.in
~/build/oofem/src/converter/converter.exec control.in      mesh.nodes mesh.voronoi
