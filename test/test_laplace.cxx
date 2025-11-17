#include <gtest/gtest.h>
#include "Laplace.hxx"
#include <gmsh.h>
#include <cmath>
#include <filesystem>

// Build a half-disk mesh (upper half of unit disk y>=0) using Gmsh and import into Laplace (which derives from Mesh)
static Laplace buildHalfDisk(double h)
{
    gmsh::initialize();
    gmsh::model::add("half_disk");
    // Geometry: circle of radius 1 centered at origin. Upper semicircle via two quarter-arc splines.
    int pL = gmsh::model::geo::addPoint(-1, 0, 0, h);
    int pT = gmsh::model::geo::addPoint(0, 1, 0, h);
    int pR = gmsh::model::geo::addPoint(1, 0, 0, h);
    int pO = gmsh::model::geo::addPoint(0, 0, 0, h);

    int arc1 = gmsh::model::geo::addCircleArc(pL, pO, pT);
    int arc2 = gmsh::model::geo::addCircleArc(pT, pO, pR);
    int base = gmsh::model::geo::addLine(pR, pL); // straight base along diameter

    int loop = gmsh::model::geo::addCurveLoop({arc1, arc2, base});
    int surf = gmsh::model::geo::addPlaneSurface({loop});

    gmsh::model::geo::synchronize();
    // Mesh size control at key points
    gmsh::vectorpair pts = { {0, pL}, {0, pT}, {0, pR} };
    gmsh::model::mesh::setSize(pts, h);
    gmsh::model::mesh::generate(2);

    Mesh M = Mesh::buildFromGmshCurrent(true);
    gmsh::finalize();

    Laplace L;
    L.nv = M.nv; L.nelems = M.nelems;
    L.verts = std::move(M.verts);
    L.elems = std::move(M.elems);
    L.vbdy = std::move(M.vbdy);
    L.edges = std::move(M.edges);
    L.elemEdges = std::move(M.elemEdges);
    L.edgeElems = std::move(M.edgeElems);
    L.areas = std::move(M.areas);
    L.vertexElems = std::move(M.vertexElems);
    L.edgeMidpoints = std::move(M.edgeMidpoints);
    return L;
}

// Boundary condition u = y on the boundary; the harmonic extension in the disk is u=y, so this is an exact solution.
static double bc_u_equals_y(double x, double y) { (void)x; return y; }

static double l2_error_y(const Laplace &L)
{
    // Compare nodal values to exact u=y at vertices.
    double err2 = 0.0, norm2 = 0.0;
    for (int i = 0; i < L.nv; ++i) {
        double ux = L.u[static_cast<size_t>(i)];
        double uy = L.verts[static_cast<size_t>(i)][1];
        double e = ux - uy; err2 += e*e; norm2 += uy*uy;
    }
    return std::sqrt(err2 / std::max(1e-16, norm2));
}

TEST(Laplace, HalfDisk_P1)
{
    Laplace L = buildHalfDisk(0.1);
    L.setDirichletFunction(&bc_u_equals_y);
    ASSERT_TRUE(L.solve(FEType::P1));
    // Write VTK and check file exists
    const std::string outP1 = "laplace_half_disk_p1.vtk";
    ASSERT_TRUE(L.writeVTK(outP1));
    ASSERT_TRUE(std::filesystem::exists(outP1));
    // The discrete solution should be close to u=y; allow modest tolerance
    double rel = l2_error_y(L);
    EXPECT_LT(rel, 5e-2);
}

TEST(Laplace, HalfDisk_CR)
{
    Laplace L = buildHalfDisk(0.1);
    L.setDirichletFunction(&bc_u_equals_y);
    ASSERT_TRUE(L.solve(FEType::CR));
    const std::string outCR = "laplace_half_disk_cr.vtk";
    ASSERT_TRUE(L.writeVTK(outCR));
    ASSERT_TRUE(std::filesystem::exists(outCR));
    double rel = l2_error_y(L);
    EXPECT_LT(rel, 7e-2);
}
