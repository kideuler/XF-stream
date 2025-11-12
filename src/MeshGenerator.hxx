#ifndef MESH_GENERATOR_HXX
#define MESH_GENERATOR_HXX

#include "RegionBuilder.hxx"
#include <string>
#include <vector>

class MeshGenerator {
public:
    // Generate a 2D mesh for given regions using Gmsh and write to mshPath.
    // h: target element size
    // Returns true on success. (Gmsh is required at build time.)
    static bool generate(const std::vector<Region>& regions,
                         const std::vector<Nurbs>& curves,
                         double h,
                         const std::string& mshPath);
};

#endif // MESH_GENERATOR_HXX
