// Implementation for Nurbs class
#include "Nurbs.hxx"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <numeric>

namespace {
// Helper to check non-decreasing vector
template <class T>
bool is_non_decreasing(const std::vector<T>& v) {
	for (std::size_t i = 1; i < v.size(); ++i) {
		if (v[i] < v[i - 1]) return false;
	}
	return true;
}

// Local helpers for interpolation (placed in anonymous namespace)
static int interp_findSpan(double u, int p, int n, const std::vector<double>& U) {
	if (u >= U[static_cast<std::size_t>(n + 1)]) return n;
	if (u <= U[static_cast<std::size_t>(p)]) return p;
	int low = p;
	int high = n + 1;
	int mid = (low + high) / 2;
	while (u < U[static_cast<std::size_t>(mid)] || u >= U[static_cast<std::size_t>(mid + 1)]) {
		if (u < U[static_cast<std::size_t>(mid)]) high = mid; else low = mid;
		mid = (low + high) / 2;
	}
	return mid;
}

static void interp_basisFuns(int span, double u, int p, const std::vector<double>& U, std::vector<double>& N) {
	N.assign(static_cast<std::size_t>(p + 1), 0.0);
	std::vector<double> left(static_cast<std::size_t>(p + 1), 0.0), right(static_cast<std::size_t>(p + 1), 0.0);
	N[0] = 1.0;
	for (int j = 1; j <= p; ++j) {
		left[static_cast<std::size_t>(j)] = u - U[static_cast<std::size_t>(span + 1 - j)];
		right[static_cast<std::size_t>(j)] = U[static_cast<std::size_t>(span + j)] - u;
		double saved = 0.0;
		for (int r = 0; r < j; ++r) {
			const double denom = right[static_cast<std::size_t>(r + 1)] + left[static_cast<std::size_t>(j - r)];
			double temp = 0.0;
			if (denom != 0.0) temp = N[static_cast<std::size_t>(r)] / denom;
			N[static_cast<std::size_t>(r)] = saved + right[static_cast<std::size_t>(r + 1)] * temp;
			saved = left[static_cast<std::size_t>(j - r)] * temp;
		}
		N[static_cast<std::size_t>(j)] = saved;
	}
}

static std::vector<double> interp_solveDense(const std::vector<double>& A_in, const std::vector<double>& bx, int Nsize) {
	std::vector<double> M = A_in; // copy
	std::vector<double> b = bx;
	for (int i = 0; i < Nsize; ++i) {
		int piv = i;
		double maxv = std::fabs(M[static_cast<std::size_t>(i * Nsize + i)]);
		for (int r = i + 1; r < Nsize; ++r) {
			double v = std::fabs(M[static_cast<std::size_t>(r * Nsize + i)]);
			if (v > maxv) { maxv = v; piv = r; }
		}
		if (maxv == 0.0) throw std::runtime_error("Interpolate: singular system");
		if (piv != i) {
			for (int c = i; c < Nsize; ++c) std::swap(M[static_cast<std::size_t>(i * Nsize + c)], M[static_cast<std::size_t>(piv * Nsize + c)]);
			std::swap(b[i], b[piv]);
		}
		const double diag = M[static_cast<std::size_t>(i * Nsize + i)];
		for (int r = i + 1; r < Nsize; ++r) {
			const double factor = M[static_cast<std::size_t>(r * Nsize + i)] / diag;
			if (factor == 0.0) continue;
			for (int c = i; c < Nsize; ++c) M[static_cast<std::size_t>(r * Nsize + c)] -= factor * M[static_cast<std::size_t>(i * Nsize + c)];
			b[r] -= factor * b[i];
		}
	}
	std::vector<double> x(static_cast<std::size_t>(Nsize), 0.0);
	for (int i = Nsize - 1; i >= 0; --i) {
		double s = b[i];
		for (int c = i + 1; c < Nsize; ++c) s -= M[static_cast<std::size_t>(i * Nsize + c)] * x[static_cast<std::size_t>(c)];
		x[static_cast<std::size_t>(i)] = s / M[static_cast<std::size_t>(i * Nsize + i)];
	}
	return x;
}
}

Nurbs::Nurbs() = default;

Nurbs::Nurbs(int degree,
			 const std::vector<Point>& controlPoints,
			 const std::vector<double>& knotVector,
			 const std::vector<double>& weights)
	: p_(degree),
	  ctrlPts_(controlPoints),
	  weights_(weights),
	  knots_(knotVector) {
	ensureWeightsSized();
}

int Nurbs::degree() const { return p_; }

int Nurbs::id() const { return id_; }

std::size_t Nurbs::dimension() const { return 2; }

std::size_t Nurbs::numControlPoints() const { return ctrlPts_.size(); }

const std::vector<Nurbs::Point>& Nurbs::controlPoints() const { return ctrlPts_; }

const std::vector<double>& Nurbs::weights() const { return weights_; }

const std::vector<double>& Nurbs::knots() const { return knots_; }

void Nurbs::setDegree(int degree) { p_ = degree; }

void Nurbs::setId(int id) { id_ = id; }

void Nurbs::setControlPoints(const std::vector<Point>& controlPoints) {
	ctrlPts_ = controlPoints;
	ensureWeightsSized();
}

void Nurbs::setWeights(const std::vector<double>& weights) {
	weights_ = weights;
	ensureWeightsSized();
}

void Nurbs::setKnots(const std::vector<double>& knotVector) { knots_ = knotVector; }

bool Nurbs::isValid(std::string* reason) const {
	// Degree non-negative
	if (p_ < 0) {
		if (reason) *reason = "Degree must be non-negative";
		return false;
	}
	// Control points non-empty and consistent dimension
	if (ctrlPts_.empty()) {
		if (reason) *reason = "Control points are empty";
		return false;
	}
	// Points are fixed-size arrays; no dimension inconsistency possible

	// Weights size
	if (weights_.size() != ctrlPts_.size()) {
		if (reason) *reason = "Weights size must match number of control points";
		return false;
	}
	for (double w : weights_) {
		if (!(w > 0.0)) { // strictly positive
			if (reason) *reason = "Weights must be strictly positive";
			return false;
		}
	}

	// Knot vector size and monotonicity
	const int n = static_cast<int>(ctrlPts_.size()) - 1;
	const int m = n + p_ + 1;
	if (static_cast<int>(knots_.size()) != m + 1) {
		if (reason) *reason = "Knot vector size must be n + p + 2";
		return false;
	}
	if (!is_non_decreasing(knots_)) {
		if (reason) *reason = "Knot vector must be non-decreasing";
		return false;
	}

	// Domain must not be degenerate
	if (uMax() <= uMin()) {
		if (reason) *reason = "Invalid parameter domain";
		return false;
	}
	return true;
}

bool Nurbs::isClamped(double tol) const {
	if (ctrlPts_.empty()) return false;
	const int n = static_cast<int>(ctrlPts_.size()) - 1;
	const int m = n + p_ + 1;
	if (static_cast<int>(knots_.size()) != m + 1) return false;

	// Check first and last p+1 knots are equal (within tol)
	const double u0 = knots_.front();
	const double um = knots_.back();
	for (int i = 0; i <= p_; ++i) {
		if (std::fabs(knots_[i] - u0) > tol) return false;
		if (std::fabs(knots_[m - i] - um) > tol) return false;
	}
	return true;
}

double Nurbs::uMin() const {
	if (knots_.empty()) return 0.0;
	return knots_.front();
}

double Nurbs::uMax() const {
	if (knots_.empty()) return 0.0;
	return knots_.back();
}

Nurbs::Point Nurbs::evaluate(double u) const {
	std::string why;
	if (!isValid(&why)) {
		throw std::runtime_error(std::string("Invalid NURBS: ") + why);
	}
	if (u < uMin() || u > uMax()) {
		throw std::runtime_error("Parameter u out of range");
	}
	const int n = static_cast<int>(ctrlPts_.size()) - 1;
	const int span = (u == uMax()) ? n : findSpan(u);

	std::vector<double> N(p_ + 1, 0.0);
	basisFuns(span, u, N);

	// Compute weighted combination
	Point C{}; // value-initialized to {0.0, 0.0}
	double wsum = 0.0;
	for (int j = 0; j <= p_; ++j) {
		const int idx = span - p_ + j;
		const double Nj = N[j];
		const double wj = weights_[static_cast<std::size_t>(idx)];
		const double coeff = Nj * wj;
		C[0] += coeff * ctrlPts_[static_cast<std::size_t>(idx)][0];
		C[1] += coeff * ctrlPts_[static_cast<std::size_t>(idx)][1];
		wsum += coeff;
	}
	if (wsum == 0.0) {
		throw std::runtime_error("Zero weight sum during evaluation");
	}
	C[0] /= wsum; C[1] /= wsum;
	return C;
}

int Nurbs::findSpan(double u) const {
	const int n = static_cast<int>(ctrlPts_.size()) - 1;
	if (u >= knots_[static_cast<std::size_t>(n + 1)]) return n;
	if (u <= knots_[static_cast<std::size_t>(p_)]) return p_;

	int low = p_;
	int high = n + 1;
	int mid = (low + high) / 2;
	while (u < knots_[static_cast<std::size_t>(mid)] || u >= knots_[static_cast<std::size_t>(mid + 1)]) {
		if (u < knots_[static_cast<std::size_t>(mid)]) {
			high = mid;
		} else {
			low = mid;
		}
		mid = (low + high) / 2;
	}
	return mid;
}

void Nurbs::basisFuns(int span, double u, std::vector<double>& N) const {
	N.assign(p_ + 1, 0.0);
	std::vector<double> left(p_ + 1, 0.0), right(p_ + 1, 0.0);
	N[0] = 1.0;
	for (int j = 1; j <= p_; ++j) {
		left[static_cast<std::size_t>(j)] = u - knots_[static_cast<std::size_t>(span + 1 - j)];
		right[static_cast<std::size_t>(j)] = knots_[static_cast<std::size_t>(span + j)] - u;
		double saved = 0.0;
		for (int r = 0; r < j; ++r) {
			const double denom = right[static_cast<std::size_t>(r + 1)] + left[static_cast<std::size_t>(j - r)];
			double temp = 0.0;
			if (denom != 0.0) temp = N[static_cast<std::size_t>(r)] / denom;
			N[static_cast<std::size_t>(r)] = saved + right[static_cast<std::size_t>(r + 1)] * temp;
			saved = left[static_cast<std::size_t>(j - r)] * temp;
		}
		N[static_cast<std::size_t>(j)] = saved;
	}
}

void Nurbs::ensureWeightsSized() {
	if (weights_.size() != ctrlPts_.size()) {
		weights_.assign(ctrlPts_.size(), 1.0);
	}
}

// Static: Interpolate points with degree p using global interpolation.
// Reference: The NURBS Book (Piegl & Tiller) Algorithm A9.1 (for B-spline), here with unit weights.
Nurbs Nurbs::interpolate(const std::vector<Point>& Q, int p) {
	if (p < 1) throw std::runtime_error("Interpolate: degree must be >= 1");
	const int mQ = static_cast<int>(Q.size());
	if (mQ < p + 1) throw std::runtime_error("Interpolate: insufficient points for given degree");

	const int n = mQ - 1; // n = number of control points - 1 (here equals points-1 for interpolation)

	// 1) Chord-length parameterization
	std::vector<double> Ubar(mQ, 0.0);
	double total = 0.0;
	for (int i = 1; i < mQ; ++i) {
		double dx = Q[i][0] - Q[i-1][0];
		double dy = Q[i][1] - Q[i-1][1];
		total += std::sqrt(dx*dx + dy*dy);
		Ubar[i] = total;
	}
	if (total > 0.0) {
		for (int i = 1; i < mQ; ++i) Ubar[i] /= total;
	} else {
		// All points coincide; return a degenerate curve
		std::vector<Point> ctrl(1, Q.front());
		std::vector<double> knots(2*p + 2, 0.0);
		for (int i = p+1; i < 2*p+2; ++i) knots[static_cast<std::size_t>(i)] = 1.0;
		return Nurbs(p, ctrl, knots);
	}

	// 2) Construct clamped knot vector by averaging (The NURBS Book eq. 9.8)
	std::vector<double> U(n + p + 2, 0.0);
	for (int j = 0; j <= p; ++j) U[static_cast<std::size_t>(j)] = 0.0;
	for (int j = n + 1; j <= n + p + 1; ++j) U[static_cast<std::size_t>(j)] = 1.0;
	for (int j = 1; j <= n - p; ++j) {
		double sum = 0.0;
		for (int i = j; i < j + p; ++i) sum += Ubar[i];
		U[static_cast<std::size_t>(j + p)] = sum / p;
	}

	// 3) Solve for control points P such that N_{i,p}(Ubar_k) P_i = Q_k
	// Build coefficient matrix A (size (n+1)x(n+1)) banded width (p+1)
	const int Nsize = n + 1;
	std::vector<double> A(static_cast<std::size_t>(Nsize * Nsize), 0.0);
	// We’ll compute spans using the same knot vector; use helper interp_findSpan

	// Fill A
	std::vector<double> Np(static_cast<std::size_t>(p + 1), 0.0);
	for (int k = 0; k <= n; ++k) {
		int span = interp_findSpan(Ubar[k], p, n, U);
		interp_basisFuns(span, Ubar[k], p, U, Np);
		for (int j = 0; j <= p; ++j) {
			int col = span - p + j;
			A[static_cast<std::size_t>(k * Nsize + col)] = Np[static_cast<std::size_t>(j)];
		}
	}


	std::vector<double> bx(static_cast<std::size_t>(Nsize), 0.0), by(static_cast<std::size_t>(Nsize), 0.0);
	for (int i = 0; i <= n; ++i) { bx[static_cast<std::size_t>(i)] = Q[static_cast<std::size_t>(i)][0]; by[static_cast<std::size_t>(i)] = Q[static_cast<std::size_t>(i)][1]; }
	std::vector<double> Px = interp_solveDense(A, bx, Nsize);
	std::vector<double> Py = interp_solveDense(A, by, Nsize);

	std::vector<Point> Pts(static_cast<std::size_t>(Nsize));
	for (int i = 0; i <= n; ++i) Pts[static_cast<std::size_t>(i)] = Point{Px[static_cast<std::size_t>(i)], Py[static_cast<std::size_t>(i)]};
	std::vector<double> W(static_cast<std::size_t>(Nsize), 1.0);

	Nurbs curve(p, Pts, U, W);
	std::string why;
	if (!curve.isValid(&why)) throw std::runtime_error(std::string("Interpolate: produced invalid curve: ") + why);
	return curve;
}

