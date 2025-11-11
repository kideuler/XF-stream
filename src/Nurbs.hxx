#ifndef __NURBS_HXX__
#define __NURBS_HXX__


#include <cstddef>
#include <string>
#include <vector>
#include <array>

// A lightweight NURBS curve container with evaluation support.
// Notes:
// - Represents a 1D parametric NURBS curve C(u) in R^dim.
// - Control points are stored as vectors of size == dimension().
// - If weights are omitted, unit weights are assumed.
// - Knot vector is non-decreasing with size m+1 where m = n + p + 1.
// - All member functions are defined in Nurbs.cxx per request.
class Nurbs {
public:
	using Point = std::array<double, 2>; // 2D point (x, y)

	// Constructors
	Nurbs();
	Nurbs(int degree,
		  const std::vector<Point>& controlPoints,
		  const std::vector<double>& knotVector,
		  const std::vector<double>& weights = {});

	// Basic metadata
	int degree() const;
	std::size_t dimension() const; // Always 2
	std::size_t numControlPoints() const;

	// Accessors
	const std::vector<Point>& controlPoints() const;
	const std::vector<double>& weights() const;
	const std::vector<double>& knots() const;

	// Mutators (re-validate on change)
	void setDegree(int degree);
	void setControlPoints(const std::vector<Point>& controlPoints);
	void setWeights(const std::vector<double>& weights);
	void setKnots(const std::vector<double>& knotVector);

	// Validation and properties
	bool isValid(std::string* reason = nullptr) const;
	bool isClamped(double tol = 0.0) const; // multiplicity p+1 at ends

	// Parameter domain [u_min, u_max] where evaluation is defined
	double uMin() const;
	double uMax() const;

	// Evaluate curve point at parameter u using Cox–de Boor
	// Throws std::runtime_error if invalid configuration or u out of range
	Point evaluate(double u) const;

	// Construct a NURBS curve of given degree that interpolates the provided points.
	// Uses a global interpolation algorithm (parameterization + knot averaging) with unit weights.
	// Throws std::runtime_error on invalid input (degree < 1, insufficient points, etc.).
	static Nurbs interpolate(const std::vector<Point>& points, int degree);

private:
	// Internal helpers
	int findSpan(double u) const; // returns span index in [p, n]
	void basisFuns(int span, double u, std::vector<double>& N) const; // size p+1
	void ensureWeightsSized();

	int p_{0};
	std::vector<Point> ctrlPts_;
	std::vector<double> weights_; // size == ctrlPts_.size()
	std::vector<double> knots_;   // non-decreasing, size m+1 where m = n+p+1
};



#endif // __NURBS_HXX__