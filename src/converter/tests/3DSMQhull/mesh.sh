~/build/oofem/src/generator/generator.exec mesh.in
mv nodes.dat mesh.nodes
qvoronoi p Fv < mesh.nodes > mesh.voronoi
~/build/oofem/src/converter/converter.exec control.in mesh.nodes mesh.voronoi
