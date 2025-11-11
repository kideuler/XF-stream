#include <gtest/gtest.h>
#include "Nurbs.hxx"
#include <cmath>

TEST(Nurbs, DefaultInvalid) {
    Nurbs c;
    std::string why;
    EXPECT_FALSE(c.isValid(&why));
}

TEST(Nurbs, ValidConstructionAndEvaluate) {
    int p = 2;
    std::vector<Nurbs::Point> P = {
        Nurbs::Point{0.0, 0.0},
        Nurbs::Point{1.0, 2.0},
        Nurbs::Point{3.0, 3.0},
        Nurbs::Point{4.0, 0.0}
    };
    std::vector<double> U = {0,0,0,1,2,2,2};
    std::vector<double> W = {1, 1, 1, 1};

    Nurbs c(p, P, U, W);
    std::string why;
    ASSERT_TRUE(c.isValid(&why)) << why;

    // Endpoints should equal first/last control points for clamped curve
    auto a = c.evaluate(c.uMin());
    auto b = c.evaluate(c.uMax());
    // Fixed 2D
    EXPECT_NEAR(a[0], P.front()[0], 1e-12);
    EXPECT_NEAR(a[1], P.front()[1], 1e-12);
    EXPECT_NEAR(b[0], P.back()[0], 1e-12);
    EXPECT_NEAR(b[1], P.back()[1], 1e-12);

    // Mid evaluation should be finite and within reasonable bounds
    auto mid = c.evaluate(1.0);
    // Fixed 2D
    EXPECT_TRUE(std::isfinite(mid[0]));
    EXPECT_TRUE(std::isfinite(mid[1]));
}

TEST(Nurbs, InterpolateCircle) {
    // Sample points on a unit circle quarter and interpolate with degree 3
    const int degree = 3;
    std::vector<Nurbs::Point> pts;
    const int samples = 9; // enough samples
    const double PI = 3.14159265358979323846;
    for (int i = 0; i < samples; ++i) {
        double t = (static_cast<double>(i) / (samples - 1)) * (PI / 2.0);
        pts.push_back(Nurbs::Point{std::cos(t), std::sin(t)});
    }
    Nurbs c = Nurbs::interpolate(pts, degree);
    std::string why;
    ASSERT_TRUE(c.isValid(&why)) << why;

    // Evaluate at parameter sites used for interpolation (approximated by chord-length from implementation)
    // We don't expose the parameterization, but evaluating at knots near interior averages should be fine; instead
    // we check the curve passes close to the original samples by searching for closest parameter with a small grid.
    auto radius_err = [&](const Nurbs::Point& p){ return std::fabs(std::sqrt(p[0]*p[0]+p[1]*p[1]) - 1.0); };
    double maxErr = 0.0;
    for (const auto& s : pts) {
        // Coarse search over u
        double best = 1e9; double bestU = c.uMin();
        for (int k = 0; k <= 200; ++k) {
            double u = c.uMin() + (c.uMax()-c.uMin()) * (static_cast<double>(k)/200.0);
            auto p = c.evaluate(u);
            double e2 = (p[0]-s[0])*(p[0]-s[0]) + (p[1]-s[1])*(p[1]-s[1]);
            if (e2 < best) { best = e2; bestU = u; }
        }
        auto pbest = c.evaluate(bestU);
        maxErr = std::max(maxErr, std::sqrt(best));
    }
    // Expect reasonably small fit error for cubic interpolation of quarter circle samples
    EXPECT_LT(maxErr, 5e-2);
}

TEST(Nurbs, WeightDefaultsToOnes) {
    int p = 1;
    std::vector<Nurbs::Point> P = { Nurbs::Point{0.0, 0.0}, Nurbs::Point{1.0, 0.0} };
    std::vector<double> U = {0,0,1,1};

    Nurbs c(p, P, U); // no weights
    std::string why;
    ASSERT_TRUE(c.isValid(&why)) << why;

    auto v = c.evaluate(0.5);
    EXPECT_NEAR(v[0], 0.5, 1e-12);
    EXPECT_NEAR(v[1], 0.0, 1e-12);
}

TEST(Nurbs, InvalidWeightsOrKnots) {
    int p = 2;
    std::vector<Nurbs::Point> P = { Nurbs::Point{0.0, 0.0}, Nurbs::Point{1.0, 0.0}, Nurbs::Point{2.0, 0.0} };
    std::vector<double> Ubad = {0,0,0,1,1,1,1};
    std::vector<double> Wbad = {1.0, 0.0, 1.0};

    Nurbs c1(p, P, Ubad);
    std::string why;
    EXPECT_FALSE(c1.isValid(&why));

    std::vector<double> U = {0,0,0,1,1,1};
    Nurbs c2(p, P, U, Wbad);
    EXPECT_FALSE(c2.isValid(&why));
}
