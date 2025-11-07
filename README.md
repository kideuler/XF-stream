This is hopefully my first successful attempt at doing Cross-field mesh generation which I feel is one of the most gorgeous and geometric pleasing ways to do quadrilateral meshing in 2D. The goal of this project is to go from some collection of NURBS geometry (with no topology) and from there get back the partitioning of the geometry into 4-sided pieces, the interface will likely be a binary with file input and output. 


The general steps of this project will be the following:

- write a small, lightweight nurbs library which takes in, control points, knot vector, degree, to construct a nurbs curve and can do interpolation.
- From a "soup" of nurbs curves construct sets of closed loops and regions.
- Transfer this geometrtic information and mesh using GMSH.
- Grab the mesh and solve the quasilinear or use the MBO method to generate the crossfields using MFEM.
- grab the data and detect and classify singularities.


Other Tasks:
- Cmake
- CI
- good command line parser and good `-h` 
- spack package so it can be installed using spack.