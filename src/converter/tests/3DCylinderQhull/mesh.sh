~/build/oofem/src/generator/generator.exec mesh.in >stdMeshGenerator.out
qvoronoi p Fv <nodes.dat > voronoi.dat;
~/build/oofem/src/converter/converter.exec control.in nodes.dat voronoi.dat >stdMeshConverter.out 
