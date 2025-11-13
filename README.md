This is hopefully my first successful attempt at doing Cross-field mesh generation which I feel is one of the most gorgeous and geometric pleasing ways to do quadrilateral meshing in 2D. The goal of this project is to go from some collection of NURBS geometry (with no topology) and from there get back the partitioning of the geometry into 4-sided pieces, the interface will likely be a binary with file input and output. 


The general steps of this project will be the following:

- write a small, lightweight nurbs library which takes in, control points, knot vector, degree, to construct a nurbs curve and can do interpolation.
- From a "soup" of nurbs curves construct sets of closed loops and regions.
- Transfer this geometrtic information and mesh using GMSH.
- Grab the mesh and solve the quasilinear or use the MBO method to generate the crossfields using MFEM.
- grab the data and detect and classify singularities.


Other Tasks:
- Cmake
- maybe some spack CI? and cacheing the build so it doesnt do everything every CI run.
- CI
- good command line parser and good `-h` 
- spack package so it can be installed using spack.


Current list of deps:
- googletest
- gmsh

Prompts:
#codebase okay, now for something also hard. I want a new hxx, cxx file which stores the mesh created from gmsh, it should be a class with all public elements and has the following:

int nelems; // number of elements.
int nv; // number of verts.
std::vector<std::array<int,3>> elems; // topology of tris.
std::vector<std::array<diuble,2>> verts; // point locations.
std::vector<int> vbdy; // vector of points indices that are on the boundary. (includes )