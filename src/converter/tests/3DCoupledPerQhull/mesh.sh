~/build/oofem/src/generator/generator.exec mesh.in
mv nodes.dat mesh.nodes
qvoronoi p Fv < mesh.nodes > mesh.voronoi

# The registered ctest runs only the SM subproblem (control.in → oofem.in,
# diffed against reference.in). The other two templates (control.tm.in and
# control.smtm.in) are kept as examples showing how the full coupled
# staggered-SMTM pipeline assembles. Once the test framework supports
# multi-target output per test directory, all three should become active.
# The SM invocation is run LAST so `oofem.in` ends up holding the SM output.
~/build/oofem/src/converter/converter.exec control.tm.in   mesh.nodes mesh.voronoi && mv oofem.in oofem.tm.in
~/build/oofem/src/converter/converter.exec control.smtm.in mesh.nodes mesh.voronoi && mv oofem.in oofem.smtm.in
~/build/oofem/src/converter/converter.exec control.in      mesh.nodes mesh.voronoi
