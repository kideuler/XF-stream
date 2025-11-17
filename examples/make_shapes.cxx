#include "Nurbs.hxx"
#include "LoopBuilder.hxx"
#include "RegionBuilder.hxx"
#include "NbsIO.hxx"
#include "NbtIO.hxx"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

static Nurbs makeLine(const Nurbs::Point& a, const Nurbs::Point& b) {
	return Nurbs(1, std::vector<Nurbs::Point>{a,b}, std::vector<double>{0,0,1,1});
}

// Exact rational quadratic quarter-circle arc around origin with radius r.
// qStart: 0 => 0->pi/2, 1 => pi/2->pi, 2 => pi->3pi/2, 3 => 3pi/2->2pi
// If reversed is true, orientation is reversed.
static Nurbs makeQuarterArcR(double r, int qStart, bool reversed = false) {
	// Cardinal endpoints
	const Nurbs::Point P[4] = { { r, 0 }, { 0, r }, { -r, 0 }, { 0, -r } };
	// Mid control for each quarter (lies on the bisector, not on the circle)
	const Nurbs::Point M[4] = { { r, r }, { -r, r }, { -r, -r }, { r, -r } };
	const int a = ((qStart % 4) + 4) % 4;
	const int b = (a + 1) % 4;
	const double w = std::sqrt(2.0) / 2.0;
	std::vector<Nurbs::Point> cps;
	std::vector<double> W;
	if (!reversed) {
		cps = { P[a], M[a], P[b] };
		W   = { 1.0, w, 1.0 };
	} else {
		cps = { P[b], M[a], P[a] };
		W   = { 1.0, w, 1.0 };
	}
	std::vector<double> U{0,0,0,1,1,1};
	return Nurbs(2, cps, U, W);
}

static void assignSequentialIds(std::vector<Nurbs>& curves) {
	for (std::size_t i = 0; i < curves.size(); ++i) curves[i].setId(static_cast<int>(i));
}

// Build an Archimedean spiral Nurbs via interpolation from radius r0 at angle a0
// to radius r1 at angle a1 (inclusive). The spiral starts at (-r0,0) when a0=pi
// and ends at (-r1,0) when a1=3*pi.
static Nurbs makeSpiral(double r0, double r1, double a0, double a1, int samples = 64, int degree = 3) {
	std::vector<Nurbs::Point> pts; pts.reserve(static_cast<std::size_t>(samples+1));
	for (int k = 0; k <= samples; ++k) {
		double t = static_cast<double>(k) / samples;
		double ang = a0 + (a1 - a0) * t;           // angle from a0 to a1
		double rad = r0 + (r1 - r0) * t;           // linearly increasing radius
		double x = rad * std::cos(ang);
		double y = rad * std::sin(ang);
		pts.push_back(Nurbs::Point{ x, y });
	}
	return Nurbs::interpolate(pts, degree);
}

int main() {
	// Shape 1: Upper half of a disk of radius 1 centered at origin
	{
		std::vector<Nurbs> curves;
	// Construct two exact rational quarter arcs: (1,0)->(0,1) and (0,1)->(-1,0)
	curves.push_back(makeQuarterArcR(1.0, 0));
	curves.push_back(makeQuarterArcR(1.0, 1));
		// Close with bottom line from (-1,0) to (1,0) to make a region (upper half-disk)
		curves.push_back(makeLine({-1.0, 0.0}, {1.0, 0.0}));

		assignSequentialIds(curves);
		// Build loops and regions
		LoopBuilder lb(1e-9);
		auto loops = lb.buildClosedLoops(curves);
		RegionBuilder rb(1e-9, 24);
		auto regions = rb.buildRegions(curves, loops);

		// Write geometry and topology
		const std::string nbs = "upper_half_disk.nbs";
		const std::string nbt = "upper_half_disk.nbt";
		std::string err;
		if (!NbsIO::writeFile(nbs, curves, &err)) { std::fprintf(stderr, "NBS write failed: %s\n", err.c_str()); return 1; }
		if (!NbtIO::writeFile(nbt, nbs, curves, loops, regions, &err)) { std::fprintf(stderr, "NBT write failed: %s\n", err.c_str()); return 1; }
	}

	// Shape 2: Snail-shell spiral region
	// Archimedean spiral from (-1,0) (angle pi, radius 1) making one full turn to (-1.5,0)
	// (angle 3*pi, radius 1.5), then a closing line from (-1.5,0) back to (-1,0).
	{
		std::vector<Nurbs> curves;
		const double pi = 3.14159265358979323846;
		// Spiral from a0=pi (start at -1,0) to a1=3*pi (end at -1.5,0)
		Nurbs spiral = makeSpiral(1.0, 1.5, pi, 3.0*pi, 64, 3);
		curves.push_back(spiral);
		// Closing segment from end (-1.5,0) to start (-1,0)
		curves.push_back(makeLine({-1.5, 0.0}, {-1.0, 0.0}));

		assignSequentialIds(curves);
		// Build loops and regions (expect a single loop region without holes)
		LoopBuilder lb(1e-9);
		auto loops = lb.buildClosedLoops(curves);
		RegionBuilder rb(1e-9, 24);
		auto regions = rb.buildRegions(curves, loops);

		const std::string nbs = "snail_shell.nbs";
		const std::string nbt = "snail_shell.nbt";
		std::string err;
		if (!NbsIO::writeFile(nbs, curves, &err)) { std::fprintf(stderr, "NBS write failed: %s\n", err.c_str()); return 1; }
		if (!NbtIO::writeFile(nbt, nbs, curves, loops, regions, &err)) { std::fprintf(stderr, "NBT write failed: %s\n", err.c_str()); return 1; }
	}

	std::printf("Wrote example shapes: upper_half_disk.{nbs,nbt}, snail_shell.{nbs,nbt}\n");
	return 0;
}

