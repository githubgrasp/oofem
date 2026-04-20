~/build/oofem/src/generator/generator.exec mesh.in
mv nodes.dat mesh.nodes
qvoronoi p Fv < mesh.nodes > mesh.voronoi

# Three separate invocations, one per subproblem; the converter writes each
# to oofem.in, which mesh.sh renames. Once the converter/test framework
# supports multi-target output, this can collapse into a single invocation.
~/build/oofem/src/converter/converter.exec control.sm.in   mesh.nodes mesh.voronoi && mv oofem.in oofem.sm.in
~/build/oofem/src/converter/converter.exec control.tm.in   mesh.nodes mesh.voronoi && mv oofem.in oofem.tm.in
~/build/oofem/src/converter/converter.exec control.smtm.in mesh.nodes mesh.voronoi && mv oofem.in oofem.smtm.in
