~/build/oofem/src/generator/generator.exec mesh.in
mv nodes.dat mesh.nodes
qdelaunay Qt i < mesh.nodes > mesh.delaunay
qvoronoi p Fv < mesh.nodes > mesh.voronoi
~/build/oofem/src/converter/converter.exec control.in mesh.nodes mesh.delaunay mesh.voronoi
