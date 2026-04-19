~/generator/generator.exe mesh.in>stdGen.out;
qdelaunay Qt i < nodes.dat >delaunay.dat
qvoronoi p Fv <nodes.dat > voronoi.dat;
~/converter/converter.exe control.in nodes.dat delaunay.dat voronoi.dat >stdMesh.out &
