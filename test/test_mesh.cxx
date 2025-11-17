#include <gtest/gtest.h>
#include "Mesh.hxx"
#include <gmsh.h>
#include <cmath>

// Helper to generate a simple square triangulation using Gmsh directly, then build Mesh.
static Mesh buildUnitSquareMesh(double h) {
    gmsh::initialize();
    gmsh::model::add("unit_square");
    // Points
    int p1 = gmsh::model::geo::addPoint(0,0,0,h);
    int p2 = gmsh::model::geo::addPoint(1,0,0,h);
    int p3 = gmsh::model::geo::addPoint(1,1,0,h);
    int p4 = gmsh::model::geo::addPoint(0,1,0,h);
    // Lines
    int l1 = gmsh::model::geo::addLine(p1,p2);
    int l2 = gmsh::model::geo::addLine(p2,p3);
    int l3 = gmsh::model::geo::addLine(p3,p4);
    int l4 = gmsh::model::geo::addLine(p4,p1);
    int cl = gmsh::model::geo::addCurveLoop({l1,l2,l3,l4});
    int surf = gmsh::model::geo::addPlaneSurface({cl});
    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(2);
    Mesh M = Mesh::buildFromGmshCurrent(true);
    gmsh::finalize();
    return M;
}

TEST(Mesh, BasicCountsAndAreas) {
    Mesh M = buildUnitSquareMesh(0.2);
    ASSERT_GT(M.nv, 0); ASSERT_GT(M.nelems, 0);
    // Total area should approximate 1.0 (unit square). Sum all element areas.
    double sumA = 0.0; for (double a : M.areas) sumA += a;
    EXPECT_NEAR(sumA, 1.0, 1e-6);
}

TEST(Mesh, EdgeAndAdjacencyConsistency) {
    Mesh M = buildUnitSquareMesh(0.25);
    // Each triangle should have 3 edge indices
    ASSERT_EQ(M.elemEdges.size(), static_cast<std::size_t>(M.nelems));
    for (const auto& ee : M.elemEdges) {
        // check edge indices are valid
        for (int e : ee) {
            ASSERT_GE(e, 0);
            ASSERT_LT(e, static_cast<int>(M.edges.size()));
        }
    }
    // edgeElems adjacency: first entry always valid element
    for (const auto& ae : M.edgeElems) {
        ASSERT_GE(ae[0], 0); ASSERT_LT(ae[0], M.nelems);
        if (ae[1] != -1) { ASSERT_GE(ae[1], 0); ASSERT_LT(ae[1], M.nelems); }
    }
}

TEST(Mesh, BoundaryVerticesNonEmpty) {
    Mesh M = buildUnitSquareMesh(0.3);
    ASSERT_FALSE(M.vbdy.empty());
    // All boundary vertices should appear in at least one boundary edge
    std::vector<char> seen(static_cast<std::size_t>(M.nv), 0);
    for (std::size_t i = 0; i < M.edgeElems.size(); ++i) {
        if (M.edgeElems[i][1] == -1) {
            auto e = M.edges[i]; seen[static_cast<std::size_t>(e[0])] = 1; seen[static_cast<std::size_t>(e[1])] = 1; }
    }
    for (int v : M.vbdy) { ASSERT_EQ(seen[static_cast<std::size_t>(v)], 1); }
}

TEST(Mesh, VertexElementAdjacency) {
    Mesh M = buildUnitSquareMesh(0.2);
    ASSERT_EQ(M.vertexElems.size(), static_cast<std::size_t>(M.nv));
    for (std::size_t v = 0; v < M.vertexElems.size(); ++v) {
        for (int ei : M.vertexElems[v]) {
            ASSERT_GE(ei, 0); ASSERT_LT(ei, M.nelems);
            auto tri = M.elems[static_cast<std::size_t>(ei)];
            ASSERT_TRUE(tri[0] == static_cast<int>(v) || tri[1] == static_cast<int>(v) || tri[2] == static_cast<int>(v));
        }
    }
}

TEST(Mesh, EdgeMidpointsInsideSegment) {
    Mesh M = buildUnitSquareMesh(0.2);
    ASSERT_EQ(M.edgeMidpoints.size(), M.edges.size());
    for (std::size_t i = 0; i < M.edges.size(); ++i) {
        auto e = M.edges[i];
        auto A = M.verts[static_cast<std::size_t>(e[0])];
        auto B = M.verts[static_cast<std::size_t>(e[1])];
        auto Mpt = M.edgeMidpoints[i];
        // midpoint property: M = (A+B)/2
        EXPECT_NEAR(Mpt[0], 0.5*(A[0]+B[0]), 1e-12);
        EXPECT_NEAR(Mpt[1], 0.5*(A[1]+B[1]), 1e-12);
    }
}
