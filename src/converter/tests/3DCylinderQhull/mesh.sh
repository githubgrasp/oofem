~/build/oofem/src/generator/generator.exec mesh.in
mv nodes.dat mesh.nodes
qvoronoi p Fv < mesh.nodes > mesh.voronoi
~/Software/genran.git/genran.exe random.in random.dat
~/build/oofem/src/converter/converter.exec control.in mesh.nodes mesh.voronoi
