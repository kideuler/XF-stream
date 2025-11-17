#include <gtest/gtest.h>
#include "NbsIO.hxx"
#include "NbtIO.hxx"
#include "LoopBuilder.hxx"
#include "RegionBuilder.hxx"
#include <cstdio>
#include <cstring>

TEST(NbsIO, WriteAndReadSimple) {
    // Build a line and a quadratic spline
    std::vector<Nurbs> curves;
    {
        std::vector<Nurbs::Point> P{ {0,0}, {1,0} };
        std::vector<double> U{0,0,1,1};
        curves.emplace_back(1, P, U);
    }
    {
        std::vector<Nurbs::Point> P{ {0,0}, {0.5,1}, {1,0} };
        std::vector<double> U{0,0,0,1,1,1};
        std::vector<double> W{1,1,1};
        curves.emplace_back(2, P, U, W);
    }
    const char* path = "io_tmp.nbs";
    std::string err;
    ASSERT_TRUE(NbsIO::writeFile(path, curves, &err)) << err;

    // Prepend a comment and ensure parser ignores it
    {
        FILE* f = std::fopen(path, "a");
        ASSERT_NE(f, nullptr);
        std::fputs("* trailing comment\n", f);
        std::fclose(f);
    }

    std::vector<Nurbs> readCurves;
    ASSERT_TRUE(NbsIO::readFile(path, readCurves, &err)) << err;
    ASSERT_EQ(readCurves.size(), curves.size());
    EXPECT_EQ(readCurves[0].degree(), 1);
    EXPECT_EQ(readCurves[1].degree(), 2);
}

TEST(NbsIO, ReadExplicitAndInterpolateForms) {
    const char* fname = "io_mix.nbs";
    FILE* f = std::fopen(fname, "wb");
    ASSERT_NE(f, nullptr);
    const char* txt =
        "* explicit quarter-circle and an interpolated cubic\n"
        "curve degree 2\n"
        "knots 0 0 0 1 1 1\n"
        "weights 1 0.7071067811865476 1\n"
        "ctrl 1 0  1 1  0 1\n"
        "endcurve\n\n"
        "# now an interpolate block\n"
        "interpolate degree 3\n"
        "points 0 0  0.3 0.1  0.6 0.8  1 1\n"
        "endcurve\n";
    std::fwrite(txt, 1, std::strlen(txt), f);
    std::fclose(f);

    std::vector<Nurbs> r;
    std::string err;
    ASSERT_TRUE(NbsIO::readFile(fname, r, &err)) << err;
    ASSERT_EQ(r.size(), 2u);
    for (auto& c : r) { std::string why; ASSERT_TRUE(c.isValid(&why)) << why; }
}

// New test for .nbt topology file referencing .nbs geometry.
// Builds a simple square (one loop, one region) and writes .nbs and .nbt then reads .nbt.
TEST(NbtIO, WriteAndReadTopology) {
    // Build square outer and a hole (two loops, one region with one inner)
    std::vector<Nurbs> curves;
    auto makeLine = [](const Nurbs::Point& a, const Nurbs::Point& b) {
        return Nurbs(1, std::vector<Nurbs::Point>{a,b}, std::vector<double>{0,0,1,1});
    };
    // Outer square
    curves.push_back(makeLine({0,0},{1,0}));
    curves.push_back(makeLine({1,0},{1,1}));
    curves.push_back(makeLine({1,1},{0,1}));
    curves.push_back(makeLine({0,1},{0,0}));
    // Inner hole
    curves.push_back(makeLine({0.25,0.25},{0.75,0.25}));
    curves.push_back(makeLine({0.75,0.25},{0.75,0.75}));
    curves.push_back(makeLine({0.75,0.75},{0.25,0.75}));
    curves.push_back(makeLine({0.25,0.75},{0.25,0.25}));

    // Write geometry file (.nbs) assigning ids automatically
    const char* nbsPath = "topo_geom.nbs";
    std::string err;
    // Assign explicit ids
    for (std::size_t i = 0; i < curves.size(); ++i) curves[i].setId(static_cast<int>(i));
    ASSERT_TRUE(NbsIO::writeFile(nbsPath, curves, &err)) << err;

    // Build loops and region
    LoopBuilder lb(1e-9); auto loops = lb.buildClosedLoops(curves);
    ASSERT_EQ(loops.size(), 2u);
    RegionBuilder rb(1e-9, 8); auto regions = rb.buildRegions(curves, loops);
    ASSERT_EQ(regions.size(), 1u);

    // Write topology file
    const char* nbtPath = "topology.nbt";
    ASSERT_TRUE(NbtIO::writeFile(nbtPath, nbsPath, curves, loops, regions, &err)) << err;

    // Read topology back
    std::string refNbs;
    std::vector<Loop> rLoops; std::vector<Region> rRegions;
    ASSERT_TRUE(NbtIO::readFile(nbtPath, curves, refNbs, rLoops, rRegions, &err)) << err;
    EXPECT_EQ(refNbs, nbsPath);
    ASSERT_EQ(rLoops.size(), loops.size());
    ASSERT_EQ(rRegions.size(), regions.size());
    // Basic structural checks
    EXPECT_FALSE(rRegions[0].outer.empty());
    EXPECT_EQ(rRegions[0].inners.size(), regions[0].inners.size());
}
