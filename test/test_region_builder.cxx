#include <gtest/gtest.h>
#include "RegionBuilder.hxx"
#include "LoopBuilder.hxx"
#include <cmath>

// Helper: create a quadratic arc from angle a0 to a1 on circle radius r
static Nurbs makeArc(double a0, double a1, double r) {
    // Use 3-point interpolation: start, mid, end
    Nurbs::Point P0{r*std::cos(a0), r*std::sin(a0)};
    Nurbs::Point P2{r*std::cos(a1), r*std::sin(a1)};
    double am = 0.5*(a0+a1);
    Nurbs::Point P1{r*std::cos(am), r*std::sin(am)};
    std::vector<Nurbs::Point> pts{P0,P1,P2};
    // Degree 2 with simple clamped knots
    int p = 2;
    std::vector<double> U{0,0,0,1,1,1};
    std::vector<double> W{1,1,1};
    return Nurbs(p, pts, U, W);
}

TEST(RegionBuilder, OuterWithHole) {
    // Four arcs for outer radius R=2 (0..pi/2, etc.)
    const double PI = 3.14159265358979323846;
    double R = 2.0;
    double r = 0.8;
    std::vector<Nurbs> curves;
    // Outer
    curves.push_back(makeArc(0.0, PI/2.0, R));
    curves.push_back(makeArc(PI/2.0, PI, R));
    curves.push_back(makeArc(PI, 1.5*PI, R));
    curves.push_back(makeArc(1.5*PI, 2.0*PI, R));
    // Inner
    curves.push_back(makeArc(0.0, PI/2.0, r));
    curves.push_back(makeArc(PI/2.0, PI, r));
    curves.push_back(makeArc(PI, 1.5*PI, r));
    curves.push_back(makeArc(1.5*PI, 2.0*PI, r));

    // Build loops from all curves (unordered) using LoopBuilder
    LoopBuilder lb(1e-5);
    auto loops = lb.buildClosedLoops(curves);
    ASSERT_EQ(loops.size(), 2u); // outer and inner

    RegionBuilder rb(1e-3, 24);
    auto regions = rb.buildRegions(curves, loops);
    ASSERT_EQ(regions.size(), 1u); // outer region with one hole
    EXPECT_EQ(regions[0].inners.size(), 1u);
    EXPECT_FALSE(regions[0].outer.empty());
}

// Helpers for axis-aligned squares using four degree-1 segments
static Nurbs makeLineSeg(const Nurbs::Point& a, const Nurbs::Point& b) {
    int p = 1; std::vector<Nurbs::Point> P{a,b}; std::vector<double> U{0,0,1,1}; return Nurbs(p,P,U);
}

static void addSquare(std::vector<Nurbs>& curves, double x0, double y0, double x1, double y1) {
    curves.push_back(makeLineSeg({x0,y0}, {x1,y0}));
    curves.push_back(makeLineSeg({x1,y0}, {x1,y1}));
    curves.push_back(makeLineSeg({x1,y1}, {x0,y1}));
    curves.push_back(makeLineSeg({x0,y1}, {x0,y0}));
}

TEST(RegionBuilder, TwoRegionsEachWithHole) {
    // Build two disjoint outer squares with inner hole squares
    std::vector<Nurbs> curves;
    // Region A outer [0,0]-[4,4], hole [1,1]-[2,2]
    addSquare(curves, 0,0, 4,4);
    addSquare(curves, 1,1, 2,2);
    // Region B outer [10,0]-[14,4], hole [11,1]-[13,2]
    addSquare(curves, 10,0, 14,4);
    addSquare(curves, 11,1, 13,2);

    LoopBuilder lb(1e-8);
    auto loops = lb.buildClosedLoops(curves);
    ASSERT_EQ(loops.size(), 4u);

    RegionBuilder rb(1e-9, 8);
    auto regions = rb.buildRegions(curves, loops);
    ASSERT_EQ(regions.size(), 2u);
    for (const auto& R : regions) {
        EXPECT_FALSE(R.outer.empty());
        EXPECT_EQ(R.inners.size(), 1u);
    }
}

TEST(RegionBuilder, HoleTouchingEdgeCountsInside) {
    // Inner square shares a vertex with outer boundary -> should be treated as inside by on-edge tolerant check
    std::vector<Nurbs> curves;
    addSquare(curves, 0,0, 4,4);      // outer
    addSquare(curves, 2,2, 3,3);      // small inner entirely inside
    addSquare(curves, 0,0, 1,1);      // this one touches outer at (0,0)

    LoopBuilder lb(1e-8);
    auto loops = lb.buildClosedLoops(curves);
    ASSERT_EQ(loops.size(), 3u);

    RegionBuilder rb(1e-8, 10);
    auto regions = rb.buildRegions(curves, loops);
    // Expect single region with two inner loops treated as holes (touching vertex counts inside)
    ASSERT_EQ(regions.size(), 1u);
    EXPECT_EQ(regions[0].inners.size(), 2u);
}

TEST(RegionBuilder, ReversedOrientationStillContains) {
    // Create an outer square normally and an inner square with reversed orientation (by swapping endpoints)
    std::vector<Nurbs> curves;
    addSquare(curves, 0,0, 5,5);
    // reversed inner: build lines reversed order
    curves.push_back(makeLineSeg({3,3}, {1,3}));
    curves.push_back(makeLineSeg({3,1}, {3,3}));
    curves.push_back(makeLineSeg({1,1}, {3,1}));
    curves.push_back(makeLineSeg({1,3}, {1,1}));

    LoopBuilder lb(1e-8);
    auto loops = lb.buildClosedLoops(curves);
    ASSERT_EQ(loops.size(), 2u);

    RegionBuilder rb(1e-9, 8);
    auto regions = rb.buildRegions(curves, loops);
    ASSERT_EQ(regions.size(), 1u);
    EXPECT_EQ(regions[0].inners.size(), 1u);
}
