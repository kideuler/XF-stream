#include <gtest/gtest.h>
#include "NbsIO.hxx"
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
