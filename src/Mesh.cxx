#include "Mesh.hxx"
#include <gmsh.h>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>

void Mesh::clear() {
    nv = 0; nelems = 0;
    verts.clear(); elems.clear(); vbdy.clear(); edges.clear(); elemEdges.clear(); edgeElems.clear(); areas.clear(); vertexElems.clear(); edgeMidpoints.clear();
}

static Mesh buildFromGmshImpl(bool includeBoundary) {
    Mesh M; M.clear();

    // Get 2D elements (triangles) from current gmsh model
    std::vector<int> elementTypes; // gmsh element type codes
    // gmsh API expects size_t for element and node tag containers
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags);

    // Node coordinates
    std::vector<std::size_t> allNodeTags; // node tags
    std::vector<double> nodeCoord;        // xyz for each tag
    std::vector<double> nodeCoordParam;   // unused parametric coords
    gmsh::model::mesh::getNodes(allNodeTags, nodeCoord, nodeCoordParam);
    std::unordered_map<int,int> nodeTagToIndex;
    M.nv = static_cast<int>(allNodeTags.size());
    M.verts.resize(M.nv);
    for (int i = 0; i < M.nv; ++i) {
        nodeTagToIndex[static_cast<int>(allNodeTags[static_cast<std::size_t>(i)])] = i;
    }
    // nodeCoord array: (x,y,z) per node in order of allNodeTags
    for (int i = 0; i < M.nv; ++i) {
        std::size_t off = static_cast<std::size_t>(3*i);
        M.verts[static_cast<std::size_t>(i)][0] = nodeCoord[off];
        M.verts[static_cast<std::size_t>(i)][1] = nodeCoord[off+1];
    }

    // Extract triangles (gmsh type 2)
    for (std::size_t t = 0; t < elementTypes.size(); ++t) {
        if (elementTypes[t] == 2) { // 3-node triangle (type 2)
            const auto& triNodes = nodeTags[t]; // flattened node tag list (size = 3 * #tris)
            int ntri = static_cast<int>(triNodes.size() / 3);
            M.elems.reserve(M.elems.size() + static_cast<std::size_t>(ntri));
            for (int k = 0; k < ntri; ++k) {
                int aTag = static_cast<int>(triNodes[static_cast<std::size_t>(3*k)]);
                int bTag = static_cast<int>(triNodes[static_cast<std::size_t>(3*k+1)]);
                int cTag = static_cast<int>(triNodes[static_cast<std::size_t>(3*k+2)]);
                int a = nodeTagToIndex[aTag];
                int b = nodeTagToIndex[bTag];
                int c = nodeTagToIndex[cTag];
                M.elems.push_back({a,b,c});
            }
        }
    }
    M.nelems = static_cast<int>(M.elems.size());

    // Compute areas
    M.areas.resize(M.elems.size());
    for (std::size_t i = 0; i < M.elems.size(); ++i) {
        auto e = M.elems[i];
    double A = Mesh::triArea(M.verts[static_cast<std::size_t>(e[0])], M.verts[static_cast<std::size_t>(e[1])], M.verts[static_cast<std::size_t>(e[2])]);
        if (A <= 0) A = -A; // ensure positive
        M.areas[i] = A;
    }

    // Build edges and connectivity
    std::unordered_map<long long,int> edgeKeyToIndex;

    M.elemEdges.resize(M.elems.size());
    for (std::size_t ei = 0; ei < M.elems.size(); ++ei) {
        auto tri = M.elems[ei];
        int v[3] = { tri[0], tri[1], tri[2] };
        for (int epos = 0; epos < 3; ++epos) {
            int a = v[epos]; int b = v[(epos+1)%3];
            long long key = Mesh::edgeKeyPair(a,b);
            auto it = edgeKeyToIndex.find(key);
            int eIdx;
            if (it == edgeKeyToIndex.end()) {
                eIdx = static_cast<int>(M.edges.size());
                M.edges.push_back(Mesh::makeEdge(a,b));
                edgeKeyToIndex.emplace(key, eIdx);
                M.edgeElems.push_back({static_cast<int>(ei), -1});
            } else {
                eIdx = it->second;
                // second element attaches
                if (M.edgeElems[static_cast<std::size_t>(eIdx)][1] != -1) {
                    // Non-manifold edge; keep but overwrite second (simplistic handling)
                    M.edgeElems[static_cast<std::size_t>(eIdx)][1] = static_cast<int>(ei);
                } else {
                    M.edgeElems[static_cast<std::size_t>(eIdx)][1] = static_cast<int>(ei);
                }
            }
            M.elemEdges[ei][epos] = eIdx;
        }
    }

    // Edge midpoints
    M.edgeMidpoints.resize(M.edges.size());
    for (std::size_t i = 0; i < M.edges.size(); ++i) {
        auto e = M.edges[i];
        M.edgeMidpoints[i] = Mesh::midpoint(M.verts[static_cast<std::size_t>(e[0])], M.verts[static_cast<std::size_t>(e[1])]);
    }

    // Boundary vertices (if requested): vertices belonging to at least one boundary edge
    if (includeBoundary) {
        std::vector<char> isBdy(static_cast<std::size_t>(M.nv), 0);
        for (std::size_t ei = 0; ei < M.edgeElems.size(); ++ei) {
            const auto& ee = M.edgeElems[ei];
            if (ee[1] == -1) { // boundary edge
                auto e = M.edges[ei];
                isBdy[static_cast<std::size_t>(e[0])] = 1;
                isBdy[static_cast<std::size_t>(e[1])] = 1;
            }
        }
        for (int i = 0; i < M.nv; ++i) if (isBdy[static_cast<std::size_t>(i)]) M.vbdy.push_back(i);
    }

    // vertexElems adjacency
    M.vertexElems.resize(static_cast<std::size_t>(M.nv));
    for (int ei = 0; ei < M.nelems; ++ei) {
        auto tri = M.elems[static_cast<std::size_t>(ei)];
        M.vertexElems[static_cast<std::size_t>(tri[0])].push_back(ei);
        M.vertexElems[static_cast<std::size_t>(tri[1])].push_back(ei);
        M.vertexElems[static_cast<std::size_t>(tri[2])].push_back(ei);
    }

    return M;
}

Mesh Mesh::buildFromGmshCurrent(bool includeBoundary) {
    return buildFromGmshImpl(includeBoundary);
}

Mesh Mesh::buildFromMshFile(const std::string& path, bool includeBoundary) {
    gmsh::initialize();
    gmsh::open(path);
    gmsh::model::mesh::removeDuplicateNodes();
    Mesh M = buildFromGmshImpl(includeBoundary);
    gmsh::finalize();
    return M;
}
