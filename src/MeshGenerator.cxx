#include "MeshGenerator.hxx"
#include "LoopBuilder.hxx"
#include <gmsh.h>

#include <stdexcept>
#include <unordered_map>
#include <cmath>

struct PointKeyHash {
    std::size_t operator()(const std::pair<long long,long long>& k) const noexcept {
        return std::hash<long long>()((k.first << 32) ^ k.second);
    }
};

static int getOrCreatePoint(double x, double y, double h, int& nextPointTag,
                            std::unordered_map<std::pair<long long,long long>, int, PointKeyHash>& map) {
    const double scale = 1e12; // tolerance bucket
    long long ix = static_cast<long long>(std::llround(x * scale));
    long long iy = static_cast<long long>(std::llround(y * scale));
    auto key = std::make_pair(ix, iy);
    auto it = map.find(key);
    if (it != map.end()) return it->second;
    int tag = nextPointTag++;
    gmsh::model::geo::addPoint(x, y, 0.0, h, tag);
    map.emplace(key, tag);
    return tag;
}

static int addCurveEntity(const Nurbs& c, bool reversed, double h,
                          int& nextPointTag, int& nextCurveTag,
                          std::unordered_map<std::pair<long long,long long>, int, PointKeyHash>& pointMap) {
    // If only 2 control points and degree 1, create a line; else create a spline through control points.
    const auto& ctrl = c.controlPoints();
    const auto& wts  = c.weights();
    std::vector<int> ptTags;
    if (!reversed) {
        for (auto& p : ctrl) ptTags.push_back(getOrCreatePoint(p[0], p[1], h, nextPointTag, pointMap));
    } else {
        for (auto it = ctrl.rbegin(); it != ctrl.rend(); ++it)
            ptTags.push_back(getOrCreatePoint((*it)[0], (*it)[1], h, nextPointTag, pointMap));
    }
    int curveTag = nextCurveTag++;
    if (ptTags.size() == 2 && c.degree() == 1) {
        gmsh::model::geo::addLine(ptTags[0], ptTags[1], curveTag);
    } else if (c.degree() == 2 && ctrl.size() == 3 && wts.size() == 3) {
        // Detect rational quadratic quarter-circle arc: weights (1, sqrt(2)/2, 1) (within tolerance)
        double targetMid = std::sqrt(2.0) / 2.0;
        double tol = 1e-6;
        bool quarterPattern = std::fabs(wts[0] - 1.0) < tol && std::fabs(wts[2] - 1.0) < tol && std::fabs(wts[1] - targetMid) < 1e-4;
        // Axis-aligned endpoints around origin: each endpoint has one coord zero, the other ±1
        struct AxisHelper {
            static bool isAxisEndpoint(const Nurbs::Point& p) {
                double a = std::fabs(p[0]);
                double b = std::fabs(p[1]);
                double tolP = 1e-6;
                return (std::fabs(a - 1.0) < tolP && b < tolP) || (std::fabs(b - 1.0) < tolP && a < tolP);
            }
        };
        bool axisEnds = AxisHelper::isAxisEndpoint(ctrl[0]) && AxisHelper::isAxisEndpoint(ctrl[2]);
        if (quarterPattern && axisEnds) {
            // Use circle arc representation. Center assumed at origin (0,0).
            int centerTag = getOrCreatePoint(0.0, 0.0, h, nextPointTag, pointMap);
            // If reversed, ptTags already reversed.
            gmsh::model::geo::addCircleArc(ptTags[0], centerTag, ptTags[ptTags.size()-1], curveTag);
        } else {
            gmsh::model::geo::addSpline(ptTags, curveTag);
        }
    } else {
        gmsh::model::geo::addSpline(ptTags, curveTag);
    }
    return curveTag;
}

static void addCurveLoopAndSurface(const std::vector<Nurbs>& curves,
                                   const Loop& L,
                                   double h,
                                   int& nextPointTag,
                                   int& nextCurveTag,
                                   int& nextCurveLoopTag,
                                   int& nextSurfaceTag,
                                   std::unordered_map<std::pair<long long,long long>, int, PointKeyHash>& pointMap,
                                   int& loopTag,
                                   int& surfTag) {
    std::vector<int> wire;
    wire.reserve(L.size());
    for (const auto& oc : L) {
        wire.push_back(addCurveEntity(curves[oc.index], oc.reversed, h, nextPointTag, nextCurveTag, pointMap));
    }
    loopTag = nextCurveLoopTag++;
    gmsh::model::geo::addCurveLoop(wire, loopTag);
    surfTag = nextSurfaceTag++;
    gmsh::model::geo::addPlaneSurface({loopTag}, surfTag);
}

bool MeshGenerator::generate(const std::vector<Region>& regions,
                             const std::vector<Nurbs>& curves,
                             double h,
                             const std::string& mshPath) {
    gmsh::initialize();
    gmsh::model::add("xfstream");
    // Tag counters (reset per generate call)
    int nextPointTag = 1;
    int nextCurveTag = 1;
    int nextCurveLoopTag = 1;
    int nextSurfaceTag = 1;
    std::vector<int> surfaceTags;
    std::unordered_map<std::pair<long long,long long>, int, PointKeyHash> pointMap;

    for (const auto& R : regions) {
        if (R.outer.empty()) continue;
        int outerLoop = 0, outerSurf = 0;
        addCurveLoopAndSurface(curves, R.outer, h,
                               nextPointTag, nextCurveTag, nextCurveLoopTag, nextSurfaceTag,
                               pointMap, outerLoop, outerSurf);
        std::vector<int> holeLoops;
        for (const auto& H : R.inners) {
            int holeLoop = 0, holeSurf = 0;
            addCurveLoopAndSurface(curves, H, h,
                                   nextPointTag, nextCurveTag, nextCurveLoopTag, nextSurfaceTag,
                                   pointMap, holeLoop, holeSurf);
            holeLoops.push_back(holeLoop);
        }
        if (!holeLoops.empty()) {
            // Remove initial surface without holes and create with holes
            gmsh::model::geo::remove({{2, outerSurf}});
            int surfWithHoles = nextSurfaceTag++;
            std::vector<int> loopsAll; loopsAll.push_back(outerLoop);
            for (int hl : holeLoops) loopsAll.push_back(hl);
            gmsh::model::geo::addPlaneSurface(loopsAll, surfWithHoles);
            surfaceTags.push_back(surfWithHoles);
        } else {
            surfaceTags.push_back(outerSurf);
        }
    }
    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(2);
    gmsh::write(mshPath);
    gmsh::finalize();
    return true;
}
