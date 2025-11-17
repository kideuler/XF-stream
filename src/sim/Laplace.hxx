#ifndef XFSTREAM_LAPLACE_HXX
#define XFSTREAM_LAPLACE_HXX

#include "Mesh.hxx"
#include <string>
#include <vector>

enum class FEType { P1, CR };

// Laplace solver on the Mesh with Dirichlet boundary conditions.
// If MFEM is available, P1 solve can use MFEM; otherwise an internal CG-based P1 solver is used.
// CR path currently assembles solution by projecting to/from P1 (placeholder for full CR assembly).
class Laplace : public Mesh {
public:
    // Solution at vertices (for CR, this holds transferred P1 values for output)
    std::vector<double> u;

    // Boundary condition: either constant or function
    bool use_const_bc = true;
    double bc_const = 0.0;
    double (*bc_func)(double, double) = nullptr; // must be static function pointer (no lambdas)

    FEType method = FEType::P1;

    Laplace();

    void setDirichletConstant(double c);
    void setDirichletFunction(double (*f)(double, double));

    // Solve -Delta u = 0 with Dirichlet boundary conditions on all boundary vertices
    bool solve(FEType type);

    // Write mesh + solution to legacy VTK (.vtk) UnstructuredGrid
    bool writeVTK(const std::string& path) const;

    // Static inline helper to evaluate BC
    static inline double evalBC(const Laplace* self, double x, double y) {
        if (!self) return 0.0;
        if (self->use_const_bc) return self->bc_const;
        if (self->bc_func) return self->bc_func(x, y);
        return 0.0;
    }

private:
    // Internal P1 assembly/solve if MFEM is not used
    bool solveP1_Internal();

    // Placeholder CR solve: fills u by projecting a CR edge-field to vertices via simple averaging
    bool solveCR_Internal();

    // Linear algebra helpers (static inline, no lambdas)
    static inline double dot(const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0.0; const std::size_t n = a.size(); for (std::size_t i = 0; i < n; ++i) s += a[i]*b[i]; return s;
    }
    static inline void axpy(double alpha, const std::vector<double>& x, std::vector<double>& y) {
        const std::size_t n = x.size(); for (std::size_t i = 0; i < n; ++i) y[i] += alpha * x[i];
    }
    static inline void scal(double alpha, std::vector<double>& x) {
        const std::size_t n = x.size(); for (std::size_t i = 0; i < n; ++i) x[i] *= alpha;
    }
};

#endif // XFSTREAM_LAPLACE_HXX
