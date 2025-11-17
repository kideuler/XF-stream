#ifndef XFSTREAM_MESH_HXX
#define XFSTREAM_MESH_HXX

#include <vector>
#include <array>
#include <string>
#include <cstddef>

// Simple triangle mesh container for P1 Lagrange and Crouzeix-Raviart (CR) elements.
// All members are public for direct access; helper static inline functions provided.
class Mesh {
public:
    // Counts
    int nv = 0;            // number of vertices
    int nelems = 0;        // number of triangle elements

    // Geometry and topology
    std::vector<std::array<double,2>> verts;        // vertex coordinates
    std::vector<std::array<int,3>>    elems;        // element to vertex indices (CCW)

    // Boundary vertex list (unique vertex indices on boundary in ascending order)
    std::vector<int> vbdy;

    // Derived connectivity
    std::vector<std::array<int,2>> edges;           // unique undirected edges (low, high)
    std::vector<std::array<int,3>> elemEdges;       // per element: edge indices corresponding to (v0,v1),(v1,v2),(v2,v0)
    std::vector<std::array<int,2>> edgeElems;       // per edge: adjacent elements (second = -1 if boundary)

    // Areas of elements (for integration, P1 gradient scaling)
    std::vector<double> areas;

    // Vertex -> incident elements
    std::vector<std::vector<int>> vertexElems;

    // Edge midpoints (useful for CR degrees of freedom)
    std::vector<std::array<double,2>> edgeMidpoints;

    // Build from current gmsh model (model should be synchronized and 2D mesh generated)
    // Returns true on success.
    static Mesh buildFromGmshCurrent(bool includeBoundary = true);

    // Build from a .msh file using gmsh API (loads file, extracts mesh, then finalizes gmsh)
    static Mesh buildFromMshFile(const std::string& path, bool includeBoundary = true);

    // Static inline helpers (no lambdas) -----------------------------------
    static inline double triArea(const std::array<double,2>& a,
                                 const std::array<double,2>& b,
                                 const std::array<double,2>& c) {
        return 0.5 * ((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]));
    }

    static inline std::array<int,2> makeEdge(int v0, int v1) {
        if (v0 < v1) return {v0,v1};
        return {v1,v0};
    }

    static inline bool edgeEqual(const std::array<int,2>& e, int a, int b) {
        return (e[0] == a && e[1] == b) || (e[0] == b && e[1] == a);
    }

    static inline std::array<double,2> midpoint(const std::array<double,2>& a,
                                                const std::array<double,2>& b) {
        return { 0.5*(a[0]+b[0]), 0.5*(a[1]+b[1]) };
    }

    // Key for hashing undirected edge (a,b) with a<b into 64-bit (assuming 32-bit ints)
    static inline long long edgeKeyPair(int a, int b) {
        if (a > b) std::swap(a,b);
        return (static_cast<long long>(a) << 32) | static_cast<long long>(b);
    }

    // Clear all data
    void clear();
};

#endif // XFSTREAM_MESH_HXX
