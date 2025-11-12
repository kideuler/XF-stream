#include <gtest/gtest.h>
#include "RegionBuilder.hxx"
#include "MeshGenerator.hxx"
#include <cstdio>
#include <cmath>

namespace {
static Nurbs makeLine(const Nurbs::Point& a, const Nurbs::Point& b) {
    return Nurbs(1, std::vector<Nurbs::Point>{a,b}, std::vector<double>{0,0,1,1});
}

static Nurbs makeQuarterArc(const Nurbs::Point& a, const Nurbs::Point& m, const Nurbs::Point& b) {
    const double w = std::sqrt(2.0) / 2.0;
    std::vector<Nurbs::Point> cps{a, m, b};
    std::vector<double> weights{1.0, w, 1.0};
    std::vector<double> U{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    Nurbs c(2, cps, U);
    c.setWeights(weights);
    return c;
}
}

TEST(MeshGenerator, GeneratesMsh) {
    // Build a simple square region with one hole
    std::vector<Nurbs> curves;
    // outer
    curves.push_back(makeLine({0,0}, {1,0}));
    curves.push_back(makeLine({1,0}, {1,1}));
    curves.push_back(makeLine({1,1}, {0,1}));
    curves.push_back(makeLine({0,1}, {0,0}));
    // inner
    curves.push_back(makeLine({0.25,0.25}, {0.75,0.25}));
    curves.push_back(makeLine({0.75,0.25}, {0.75,0.75}));
    curves.push_back(makeLine({0.75,0.75}, {0.25,0.75}));
    curves.push_back(makeLine({0.25,0.75}, {0.25,0.25}));

    LoopBuilder lb(1e-9);
    auto loops = lb.buildClosedLoops(curves);

    RegionBuilder rb(1e-9, 12);
    auto regions = rb.buildRegions(curves, loops);

    // Generate mesh
    const char* path = "test_out.msh";
    bool ok = MeshGenerator::generate(regions, curves, 0.1, path);
    ASSERT_TRUE(ok);
    // file existence check (basic)
    FILE* f = std::fopen(path, "rb");
    ASSERT_NE(f, nullptr);
    std::fclose(f);
}

TEST(MeshGenerator, GeneratesMshCircle) {
    // Build a unit circle (radius 1) as 4 quadratic NURBS arcs (quarter circles)
    std::vector<Nurbs> curves;
    // Define 4 quarter arcs counter-clockwise starting at (1,0)
    curves.push_back(makeQuarterArc({1,0}, {1,1}, {0,1}));
    curves.push_back(makeQuarterArc({0,1}, {-1,1}, {-1,0}));
    curves.push_back(makeQuarterArc({-1,0}, {-1,-1}, {0,-1}));
    curves.push_back(makeQuarterArc({0,-1}, {1,-1}, {1,0}));

    LoopBuilder lb(1e-9);
    auto loops = lb.buildClosedLoops(curves);
    ASSERT_EQ(loops.size(), 1u);

    RegionBuilder rb(1e-9, 32);
    auto regions = rb.buildRegions(curves, loops);
    ASSERT_EQ(regions.size(), 1u);

    const char* path = "test_circle.msh";
    bool ok = MeshGenerator::generate(regions, curves, 0.05, path);
    ASSERT_TRUE(ok);
    FILE* f = std::fopen(path, "rb");
    ASSERT_NE(f, nullptr);
    std::fclose(f);
}

