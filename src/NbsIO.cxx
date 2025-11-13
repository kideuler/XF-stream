#include "NbsIO.hxx"

#include <fstream>
#include <sstream>
#include <cctype>
#include <stdexcept>

namespace {
static inline std::string trim(const std::string& s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b-1]))) --b;
    return s.substr(a, b - a);
}

static inline bool startsWith(const std::string& s, const std::string& p) {
    return s.size() >= p.size() && std::equal(p.begin(), p.end(), s.begin());
}

static void splitTokens(const std::string& line, std::vector<std::string>& out) {
    out.clear();
    std::istringstream iss(line);
    std::string tok;
    while (iss >> tok) out.push_back(tok);
}
}

namespace NbsIO {

//
// NBS text format (2D)
// Lines starting with '*' are comments. Inline comments after '#' are ignored.
// Whitespace separated tokens.
// Curve forms:
// 1) Explicit:
//    curve degree <p>
//    knots <k0> <k1> ... <km>
//    weights <w0> <w1> ... <wn>
//    ctrl <x0> <y0>  <x1> <y1> ... <xn> <yn>
//    endcurve
//
// 2) Interpolate:
//    interpolate degree <p>
//    points <x0> <y0>  <x1> <y1> ... <xN> <yN>
//    endcurve
//

static void resetCurve(int& degree, std::vector<double>& knots, std::vector<double>& weights,
                       std::vector<Nurbs::Point>& ctrl, std::vector<Nurbs::Point>& points) {
    degree = -1; knots.clear(); weights.clear(); ctrl.clear(); points.clear();
}

static void finishExplicitOrThrow(int degree, const std::vector<double>& knots, const std::vector<double>& weights,
                                  const std::vector<Nurbs::Point>& ctrl, std::vector<Nurbs>& out) {
    if (degree < 0 || ctrl.empty() || knots.empty()) throw std::runtime_error("NBS: incomplete explicit curve");
    if (!weights.empty() && weights.size() != ctrl.size()) throw std::runtime_error("NBS: weights size mismatch");
    Nurbs c(degree, ctrl, knots, weights);
    std::string why; if (!c.isValid(&why)) throw std::runtime_error(std::string("NBS: invalid explicit curve: ")+why);
    out.push_back(std::move(c));
}

static void finishInterpOrThrow(int degree, const std::vector<Nurbs::Point>& points, std::vector<Nurbs>& out) {
    if (degree < 0 || points.size() < static_cast<std::size_t>(degree+1)) throw std::runtime_error("NBS: insufficient points for interpolation");
    Nurbs c = Nurbs::interpolate(points, degree);
    out.push_back(std::move(c));
}

bool readFile(const std::string& path, std::vector<Nurbs>& out) {
    std::ifstream ifs(path);
    if (!ifs) return false;
    std::string line;
    std::vector<std::string> toks;
    enum class State { Idle, InCurveExplicit, InCurveInterp };
    State st = State::Idle;
    int degree = -1;
    std::vector<double> knots, weights;
    std::vector<Nurbs::Point> ctrl;
    std::vector<Nurbs::Point> points;
    resetCurve(degree, knots, weights, ctrl, points);

    while (std::getline(ifs, line)) {
        // Strip inline comments after '#'
        auto hashPos = line.find('#');
        if (hashPos != std::string::npos) line = line.substr(0, hashPos);
        std::string t = trim(line);
        if (t.empty()) continue;
        if (t[0] == '*') continue; // full-line comment

        splitTokens(t, toks);
        if (toks.empty()) continue;
        if (st == State::Idle) {
            if (toks[0] == "curve") {
                st = State::InCurveExplicit; resetCurve(degree, knots, weights, ctrl, points);
                if (toks.size() >= 3 && toks[1] == "degree") degree = std::stoi(toks[2]);
            } else if (toks[0] == "interpolate") {
                st = State::InCurveInterp; resetCurve(degree, knots, weights, ctrl, points);
                if (toks.size() >= 3 && toks[1] == "degree") degree = std::stoi(toks[2]);
            } else {
                throw std::runtime_error("NBS: expected 'curve' or 'interpolate'");
            }
        } else if (st == State::InCurveExplicit) {
            if (toks[0] == "degree" && toks.size() >= 2) {
                degree = std::stoi(toks[1]);
            } else if (toks[0] == "knots") {
                for (std::size_t i = 1; i < toks.size(); ++i) knots.push_back(std::stod(toks[i]));
            } else if (toks[0] == "weights") {
                for (std::size_t i = 1; i < toks.size(); ++i) weights.push_back(std::stod(toks[i]));
            } else if (toks[0] == "ctrl") {
                if ((toks.size()-1) % 2 != 0) throw std::runtime_error("NBS: ctrl requires pairs of x y");
                for (std::size_t i = 1; i+1 < toks.size(); i += 2) ctrl.push_back(Nurbs::Point{std::stod(toks[i]), std::stod(toks[i+1])});
            } else if (toks[0] == "endcurve") {
                finishExplicitOrThrow(degree, knots, weights, ctrl, out); st = State::Idle; resetCurve(degree, knots, weights, ctrl, points);
            } else {
                throw std::runtime_error("NBS: unknown token in explicit curve");
            }
        } else if (st == State::InCurveInterp) {
            if (toks[0] == "degree" && toks.size() >= 2) {
                degree = std::stoi(toks[1]);
            } else if (toks[0] == "points") {
                if ((toks.size()-1) % 2 != 0) throw std::runtime_error("NBS: points requires pairs of x y");
                for (std::size_t i = 1; i+1 < toks.size(); i += 2) points.push_back(Nurbs::Point{std::stod(toks[i]), std::stod(toks[i+1])});
            } else if (toks[0] == "endcurve") {
                finishInterpOrThrow(degree, points, out); st = State::Idle; resetCurve(degree, knots, weights, ctrl, points);
            } else {
                throw std::runtime_error("NBS: unknown token in interpolate curve");
            }
        }
    }
    if (st != State::Idle) throw std::runtime_error("NBS: unterminated curve block");
    return true;
}

bool readFile(const std::string& path, std::vector<Nurbs>& out, std::string* errorMessage) {
    try {
        return readFile(path, out);
    } catch (const std::exception& e) {
        if (errorMessage) *errorMessage = e.what();
        return false;
    }
}

bool writeFile(const std::string& path, const std::vector<Nurbs>& curves) {
    std::ofstream ofs(path);
    if (!ofs) return false;
    ofs << "* XF-stream NBS 2D format\n";
    for (const auto& c : curves) {
        ofs << "curve degree " << c.degree() << "\n";
        ofs << "knots";
        for (double u : c.knots()) ofs << ' ' << u;
        ofs << "\nweights";
        for (double w : c.weights()) ofs << ' ' << w;
        ofs << "\nctrl";
        for (const auto& p : c.controlPoints()) ofs << ' ' << p[0] << ' ' << p[1];
        ofs << "\nendcurve\n\n";
    }
    return true;
}

bool writeFile(const std::string& path, const std::vector<Nurbs>& curves, std::string* errorMessage) {
    try {
        return writeFile(path, curves);
    } catch (const std::exception& e) {
        if (errorMessage) *errorMessage = e.what();
        return false;
    }
}

} // namespace NbsIO
