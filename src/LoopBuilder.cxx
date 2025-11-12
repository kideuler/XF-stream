#include "LoopBuilder.hxx"

#include <cmath>
#include <limits>

static inline double sqr(double x) { return x * x; }

static inline double dist2(const Nurbs::Point& a, const Nurbs::Point& b) {
    return sqr(a[0] - b[0]) + sqr(a[1] - b[1]);
}

namespace {
static bool findNextCurve(const std::vector<Nurbs>& curves,
                          const std::vector<char>& used,
                          double tol,
                          const Nurbs::Point& tail,
                          std::size_t& choiceIdx,
                          bool& choiceRev) {
    const std::size_t n = curves.size();
    double bestD2 = std::numeric_limits<double>::infinity();
    std::size_t bestIdx = n;
    bool bestRev = false;
    for (std::size_t i = 0; i < n; ++i) {
        if (used[i]) continue;
        const auto s = curves[i].evaluate(curves[i].uMin());
        const auto e = curves[i].evaluate(curves[i].uMax());
        double d2s = dist2(tail, s);
        double d2e = dist2(tail, e);
        if (d2s < bestD2 && d2s <= tol * tol) { bestD2 = d2s; bestIdx = i; bestRev = false; }
        if (d2e < bestD2 && d2e <= tol * tol) { bestD2 = d2e; bestIdx = i; bestRev = true; }
    }
    if (bestIdx == n) return false;
    choiceIdx = bestIdx; choiceRev = bestRev; return true;
}
} // anonymous namespace



Nurbs::Point LoopBuilder::startPoint(const Nurbs& c, bool reversed) {
    return reversed ? c.evaluate(c.uMax()) : c.evaluate(c.uMin());
}

Nurbs::Point LoopBuilder::endPoint(const Nurbs& c, bool reversed) {
    return reversed ? c.evaluate(c.uMin()) : c.evaluate(c.uMax());
}

bool LoopBuilder::near(const Nurbs::Point& a, const Nurbs::Point& b, double tol) {
    return dist2(a, b) <= tol * tol;
}

std::vector<Loop> LoopBuilder::buildClosedLoops(const std::vector<Nurbs>& curves,
                                                std::vector<std::size_t>* leftovers) const {
    const std::size_t n = curves.size();
    std::vector<char> used(n, 0);
    std::vector<Loop> loops;


    for (std::size_t i = 0; i < n; ++i) {
        if (used[i]) continue;
        // Start a new loop from curve i; try both orientations greedily
        for (int orient = 0; orient < 2; ++orient) {
            if (used[i]) break; // might be consumed by previous attempt
            Loop current;
            current.push_back({i, orient == 1});
            used[i] = 1;

            const auto head = startPoint(curves[i], current.front().reversed);
            auto tail = endPoint(curves[i], current.front().reversed);

            // Greedily append until we either close the loop or get stuck
            bool closed = false;
            while (true) {
                if (near(tail, head, tol_)) { closed = true; break; }
                std::size_t nextIdx; bool nextRev;
                if (!findNextCurve(curves, used, tol_, tail, nextIdx, nextRev)) break;
                current.push_back({nextIdx, nextRev});
                used[nextIdx] = 1;
                tail = endPoint(curves[nextIdx], nextRev);
            }

            if (closed) {
                loops.push_back(std::move(current));
                break; // done with this seed
            } else {
                // Rollback: mark used in this attempt as unused
                for (const auto& oc : current) used[oc.index] = 0;
            }
        }
    }

    if (leftovers) {
        leftovers->clear();
        for (std::size_t i = 0; i < n; ++i) if (!used[i]) leftovers->push_back(i);
    }
    return loops;
}
