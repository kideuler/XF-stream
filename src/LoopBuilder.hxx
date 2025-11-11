#ifndef LOOP_BUILDER_HXX
#define LOOP_BUILDER_HXX

#include "Nurbs.hxx"
#include <vector>
#include <utility>
#include <optional>
#include <cstddef>

// LoopBuilder: stitches unordered NURBS curves into closed loops.
// A loop is a sequence of oriented curve indices (index into input vector, and a flag
// indicating whether the curve is reversed) such that the end of each segment
// meets the start of the next (within a tolerance) and final end meets first start.
//
// For simplicity, we treat a curve's start point as evaluate(uMin()) and end point as evaluate(uMax()).
// Reversing swaps these. We assume curves are valid and 2D.
struct OrientedCurve {
    std::size_t index; // original curve index
    bool reversed;     // true if used in reversed orientation
};

using Loop = std::vector<OrientedCurve>;

class LoopBuilder {
public:
    static Nurbs::Point startPoint(const Nurbs& c, bool reversed);
    static Nurbs::Point endPoint(const Nurbs& c, bool reversed);

    explicit LoopBuilder(double tol = 1e-8) : tol_(tol) {}

    // Build maximal set of closed loops. Each input curve used at most once.
    // Curves that cannot be stitched are ignored (returned in leftovers indices if requested).
    std::vector<Loop> buildClosedLoops(const std::vector<Nurbs>& curves,
                                       std::vector<std::size_t>* leftovers = nullptr) const;

private:
    double tol_;
    static bool near(const Nurbs::Point& a, const Nurbs::Point& b, double tol);
};

#endif // LOOP_BUILDER_HXX
