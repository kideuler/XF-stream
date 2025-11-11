#ifndef REGION_BUILDER_HXX
#define REGION_BUILDER_HXX

#include "LoopBuilder.hxx" // for Loop and OrientedCurve
#include "Nurbs.hxx"
#include <vector>

struct Region {
    // Outer loop of the region; if empty, the region has no outer (rare)
    Loop outer;
    // Inner loops (holes)
    std::vector<Loop> inners;
};

class RegionBuilder {
public:
    explicit RegionBuilder(double tol = 1e-8, int samplesPerCurve = 16)
        : tol_(tol), samplesPerCurve_(samplesPerCurve) {}

    // Build regions from loops; each loop is defined over the curves array (by indices in OrientedCurve)
    std::vector<Region> buildRegions(const std::vector<Nurbs>& curves,
                                     const std::vector<Loop>& loops) const;

private:
    double tol_;
    int samplesPerCurve_;

    // Sample points along a loop to form a polygon approximation
    std::vector<Nurbs::Point> sampleLoop(const std::vector<Nurbs>& curves, const Loop& L) const;

    // Geometric predicates
    static bool pointInPolygon(const std::vector<Nurbs::Point>& poly, const Nurbs::Point& p, double tol);
    static bool onSegmentTol(const Nurbs::Point& a, const Nurbs::Point& b, const Nurbs::Point& p, double tol);
};

#endif // REGION_BUILDER_HXX
