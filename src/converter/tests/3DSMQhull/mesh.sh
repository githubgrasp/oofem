~/build/oofem/src/generator/generator.exec mesh.in;
qvoronoi p Fv <mesh.nodes > mesh.voronoi;
~/build/oofem/src/converter/converter.exec control.in mesh.nodes mesh.voronoi;
