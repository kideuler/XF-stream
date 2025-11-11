#include "RegionBuilder.hxx"

#include <cmath>
#include <algorithm>

// Sample each curve segment with 'samplesPerCurve_' points (uniform in parameter range).
std::vector<Nurbs::Point> RegionBuilder::sampleLoop(const std::vector<Nurbs>& curves, const Loop& L) const {
    std::vector<Nurbs::Point> pts;
    if (L.empty()) return pts;
    for (std::size_t idxSeg = 0; idxSeg < L.size(); ++idxSeg) {
        const auto& oc = L[idxSeg];
        const auto& c = curves[oc.index];
        double u0 = oc.reversed ? c.uMax() : c.uMin();
        double u1 = oc.reversed ? c.uMin() : c.uMax();
        for (int s = 0; s < samplesPerCurve_; ++s) {
            // Avoid duplicating the end point of this segment if it's the start of the next
            if (s == samplesPerCurve_ - 1 && idxSeg + 1 < L.size()) continue;
            double t = (samplesPerCurve_ == 1) ? 0.0 : static_cast<double>(s) / (samplesPerCurve_ - 1);
            double u = u0 + (u1 - u0) * t;
            pts.push_back(c.evaluate(u));
        }
    }
    return pts;
}

bool RegionBuilder::onSegmentTol(const Nurbs::Point& a, const Nurbs::Point& b, const Nurbs::Point& p, double tol) {
    double vx = b[0] - a[0];
    double vy = b[1] - a[1];
    double wx = p[0] - a[0];
    double wy = p[1] - a[1];
    double c1 = vx * wx + vy * wy;
    if (c1 < -tol) return false;
    double c2 = vx * vx + vy * vy;
    if (c1 - c2 > tol) return false;
    double cross = std::fabs(vx * wy - vy * wx);
    return cross <= tol * std::sqrt(c2);
}

bool RegionBuilder::pointInPolygon(const std::vector<Nurbs::Point>& poly, const Nurbs::Point& p, double tol) {
    if (poly.size() < 3) return false;
    // Ray casting with tolerance on edges
    int crossings = 0;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        const auto& a = poly[i];
        const auto& b = poly[(i + 1) % poly.size()];
        if (onSegmentTol(a, b, p, tol)) return true; // treat on-edge as inside
        bool condY = ((a[1] <= p[1]) && (b[1] > p[1])) || ((a[1] > p[1]) && (b[1] <= p[1]));
        if (condY) {
            double xInt = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1] + 1e-20);
            if (xInt > p[0]) crossings++;
        }
    }
    return (crossings % 2) == 1;
}

std::vector<Region> RegionBuilder::buildRegions(const std::vector<Nurbs>& curves, const std::vector<Loop>& loops) const {
    // Sample all loops
    struct LoopGeom { Loop loop; std::vector<Nurbs::Point> poly; Nurbs::Point repr; };
    std::vector<LoopGeom> geoms;
    geoms.reserve(loops.size());
    for (const auto& L : loops) {
        auto poly = sampleLoop(curves, L);
        if (poly.empty()) continue;
        // Representative point: polygon centroid (robust for convex shapes; adequate for our tests)
        double cx = 0.0, cy = 0.0;
        for (const auto& p : poly) { cx += p[0]; cy += p[1]; }
        cx /= static_cast<double>(poly.size());
        cy /= static_cast<double>(poly.size());
        geoms.push_back(LoopGeom{L, poly, Nurbs::Point{cx, cy}});
    }

    // Determine containment graph: loop A contains loop B if B's representative point lies inside A.
    std::size_t m = geoms.size();
    auto polyArea = [&](const std::vector<Nurbs::Point>& poly){
        double a = 0.0; for (std::size_t k = 0; k < poly.size(); ++k){ const auto& A=poly[k]; const auto& B=poly[(k+1)%poly.size()]; a += A[0]*B[1]-B[0]*A[1]; } return std::fabs(a)*0.5; };
    std::vector<double> areas(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) areas[i] = polyArea(geoms[i].poly);

    struct Box { double minx, miny, maxx, maxy; };
    auto computeBox = [](const std::vector<Nurbs::Point>& poly){
        Box b{poly[0][0], poly[0][1], poly[0][0], poly[0][1]};
        for (const auto& p : poly) { b.minx = std::min(b.minx, p[0]); b.miny = std::min(b.miny, p[1]); b.maxx = std::max(b.maxx, p[0]); b.maxy = std::max(b.maxy, p[1]); }
        return b;
    };
    std::vector<Box> boxes(m);
    for (std::size_t i = 0; i < m; ++i) boxes[i] = computeBox(geoms[i].poly);
    auto boxContains = [&](const Box& A, const Box& B){
        return A.minx - tol_ <= B.minx && A.miny - tol_ <= B.miny && A.maxx + tol_ >= B.maxx && A.maxy + tol_ >= B.maxy;
    };
    std::vector<std::vector<std::size_t>> children(m); // indices of loops contained by i
    std::vector<int> parent(m, -1);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) if (i != j) {
            // Only allow larger area loops to contain smaller ones to avoid mutual containment cycles
            if (areas[i] > areas[j] && boxContains(boxes[i], boxes[j]) && pointInPolygon(geoms[i].poly, geoms[j].repr, tol_)) {
                // i contains j; assign parent if more outer than existing
                if (parent[j] == -1) parent[j] = static_cast<int>(i);
                else {
                    double areaOld = areas[static_cast<std::size_t>(parent[j])];
                    double areaNew = areas[i];
                    if (areaNew < areaOld) parent[j] = static_cast<int>(i);
                }
            }
        }
    }

    // Build regions: a loop with parent == -1 is an outer; loops whose parent is that outer become holes.
    std::vector<Region> regions;
    for (std::size_t i = 0; i < m; ++i) {
        if (parent[i] == -1) {
            Region R; R.outer = geoms[i].loop;
            // Collect direct children (holes). We could also filter nested holes but keep simple: only direct ones.
            for (std::size_t j = 0; j < m; ++j) if (parent[j] == static_cast<int>(i)) R.inners.push_back(geoms[j].loop);
            regions.push_back(std::move(R));
        }
    }
    // If no explicit outers found (due to tolerance or degenerate input), promote largest loops as outers by area and attach smaller ones as holes if they were parentless
    if (regions.empty() && m > 0) {
        // Build regions by grouping loops that don't contain each other's centroids.
        std::vector<char> assigned(m, 0);
        for (std::size_t i = 0; i < m; ++i) {
            if (assigned[i]) continue;
            // Start a region with loop i as outer
            Region R; R.outer = geoms[i].loop; assigned[i] = 1;
            // Any loops contained by this (by centroid) become holes; disjoint loops will start their own region
            for (std::size_t j = 0; j < m; ++j) {
                if (assigned[j] || j == i) continue;
                if (areas[i] > areas[j] && boxContains(boxes[i], boxes[j]) && pointInPolygon(geoms[i].poly, geoms[j].repr, tol_)) {
                    R.inners.push_back(geoms[j].loop); assigned[j] = 1;
                }
            }
            regions.push_back(std::move(R));
        }
    }
    // Sort regions by decreasing area for consistency
    std::sort(regions.begin(), regions.end(), [&](const Region& a, const Region& b){
        // approximate area of outer polygons
        auto areaOf = [&](const Loop& L){
            auto poly = sampleLoop(curves, L);
            double A = 0.0; for (std::size_t k = 0; k < poly.size(); ++k){ const auto& P=poly[k]; const auto& Q=poly[(k+1)%poly.size()]; A += P[0]*Q[1]-Q[0]*P[1]; } return std::fabs(A)*0.5; };
        return areaOf(a.outer) > areaOf(b.outer);
    });
    return regions;
}
