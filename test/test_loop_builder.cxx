#include <gtest/gtest.h>
#include "LoopBuilder.hxx"
#include <cmath>

static Nurbs makeLine(const Nurbs::Point& a, const Nurbs::Point& b) {
    // Degree 1 line with two control points and clamped knots {0,0,1,1}
    int p = 1;
    std::vector<Nurbs::Point> P = { a, b };
    std::vector<double> U = {0.0, 0.0, 1.0, 1.0};
    return Nurbs(p, P, U);
}

TEST(LoopBuilder, BuildsSingleTriangle) {
    auto c1 = makeLine({0.0, 0.0}, {1.0, 0.0});
    auto c2 = makeLine({1.0, 0.0}, {0.5, std::sqrt(3.0)/2.0});
    auto c3 = makeLine({0.5, std::sqrt(3.0)/2.0}, {0.0, 0.0});

    std::vector<Nurbs> curves{c2, c3, c1}; // intentionally unordered

    LoopBuilder lb(1e-10);
    auto loops = lb.buildClosedLoops(curves);
    ASSERT_EQ(loops.size(), 1u);
    const auto& L = loops[0];
    ASSERT_EQ(L.size(), 3u);
    // Verify closed by evaluating endpoints
    auto head = (L.front().reversed ? curves[L.front().index].evaluate(curves[L.front().index].uMax())
                                    : curves[L.front().index].evaluate(curves[L.front().index].uMin()));
    auto tail = (L.back().reversed ? curves[L.back().index].evaluate(curves[L.back().index].uMin())
                                   : curves[L.back().index].evaluate(curves[L.back().index].uMax()));
    double dx = head[0] - tail[0];
    double dy = head[1] - tail[1];
    EXPECT_NEAR(dx*dx + dy*dy, 0.0, 1e-16);
}

TEST(LoopBuilder, MultipleLoopsAndLeftovers) {
    // Two squares and one stray segment
    auto s1a = makeLine({0,0}, {1,0});
    auto s1b = makeLine({1,0}, {1,1});
    auto s1c = makeLine({1,1}, {0,1});
    auto s1d = makeLine({0,1}, {0,0});

    auto s2a = makeLine({2,0}, {3,0});
    auto s2b = makeLine({3,0}, {3,1});
    auto s2c = makeLine({3,1}, {2,1});
    auto s2d = makeLine({2,1}, {2,0});

    auto stray = makeLine({10,10}, {11,10});

    std::vector<Nurbs> curves{ s1a,s2d,s1c,stray,s2b,s2a,s1d,s1b,s2c };
    LoopBuilder lb(1e-12);
    std::vector<std::size_t> leftovers;
    auto loops = lb.buildClosedLoops(curves, &leftovers);

    ASSERT_EQ(loops.size(), 2u);
    EXPECT_EQ(leftovers.size(), 1u);
}
