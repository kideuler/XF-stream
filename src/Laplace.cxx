#include "Laplace.hxx"
#include <fstream>
#include <cmath>
#include <queue>
#include <mfem.hpp>

Laplace::Laplace() = default;

void Laplace::setDirichletConstant(double c) {
    use_const_bc = true; bc_const = c; bc_func = nullptr;
}
void Laplace::setDirichletFunction(double (*f)(double,double)) {
    use_const_bc = false; bc_func = f;
}

bool Laplace::solve(FEType type) {
    method = type;
    if (nv == 0 || nelems == 0) return false;
    u.assign(static_cast<std::size_t>(nv), 0.0);
    switch (type) {
        case FEType::P1: return solveP1_Internal();
        case FEType::CR: return solveCR_Internal();
        default: return false;
    }
}

// Simple P1 Laplace with Dirichlet on boundary: assemble stiffness matrix in CSR-like triplets and solve by CG with diagonal preconditioner.
// This avoids MFEM dependency for portability, but can be replaced with MFEM-based solve later.
bool Laplace::solveP1_Internal() {
    using namespace mfem;
    // Build an MFEM Mesh from our triangle mesh
    mfem::Mesh m(2, nv, nelems, /*NBdr*/0, /*spaceDim*/2);
    for (int i = 0; i < nv; ++i) { m.AddVertex(verts[static_cast<std::size_t>(i)][0], verts[static_cast<std::size_t>(i)][1]); }
    for (int e = 0; e < nelems; ++e) {
        int vi[3] = { elems[static_cast<std::size_t>(e)][0], elems[static_cast<std::size_t>(e)][1], elems[static_cast<std::size_t>(e)][2] };
        m.AddTriangle(vi, /*attr*/1);
    }
    m.FinalizeTriMesh(/*gen_edges*/1, /*refine*/0, /*fix_orientation*/true);

    // H1 P1 space
    H1_FECollection fec(1, m.Dimension());
    FiniteElementSpace fes(&m, &fec);

    // Essential boundary dofs: mark all boundary vertices
    Array<int> ess_tdof_list;
    Array<int> ess_bdr;
    if (m.bdr_attributes.Size() == 0) { m.GenerateBoundaryElements(); }
    int nb = m.bdr_attributes.Size() ? m.bdr_attributes.Max() : 1;
    ess_bdr.SetSize(nb); ess_bdr = 1; // all boundaries essential
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // Coefficient for BC
    class DirichletCoeff : public Coefficient {
        const Laplace* L;
    public: explicit DirichletCoeff(const Laplace* l): L(l) {}
        virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) {
            Vector x; T.Transform(ip, x); return Laplace::evalBC(L, x(0), x(1));
        }
    } dbc(this);

    // Setup system: (grad u, grad v) = 0
    BilinearForm a(&fes);
    a.AddDomainIntegrator(new DiffusionIntegrator);
    a.Assemble();

    GridFunction x(&fes); x = 0.0;
    // Apply boundary conditions by interpolation
    x.ProjectBdrCoefficient(dbc, ess_bdr);

    // Form linear system
    LinearForm b(&fes); b.Assemble();
    SparseMatrix A; Vector X, B;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

    // Solve with CG + Jacobi
    GSSmoother M(A);
    PCG(A, M, B, X, 1, 2000, 1e-12, 0.0);

    // Recover solution
    a.RecoverFEMSolution(X, b, x);

    // Copy to vertex order using vertex dof map
    for (int i = 0; i < nv; ++i) {
        mfem::Array<int> vdofs; fes.GetVertexDofs(i, vdofs);
        double val = 0.0; int cnt = 0;
        for (int k = 0; k < vdofs.Size(); ++k) {
            int di = vdofs[k];
            if (di < 0) di = -1 - di; // handle sign/constrained dofs
            if (di >= 0 && di < x.Size()) { val += x(di); cnt++; }
        }
        u[static_cast<std::size_t>(i)] = (cnt>0) ? (val/cnt) : 0.0;
    }
    return true;
}

// CR placeholder: produce an edge-based constant field (zero), then average to vertices for output
bool Laplace::solveCR_Internal() {
    using namespace mfem;
    mfem::Mesh m(2, nv, nelems, /*NBdr*/0, /*spaceDim*/2);
    for (int i = 0; i < nv; ++i) { m.AddVertex(verts[static_cast<std::size_t>(i)][0], verts[static_cast<std::size_t>(i)][1]); }
    for (int e = 0; e < nelems; ++e) { int vi[3] = { elems[static_cast<std::size_t>(e)][0], elems[static_cast<std::size_t>(e)][1], elems[static_cast<std::size_t>(e)][2] }; m.AddTriangle(vi, 1); }
    m.FinalizeTriMesh(1,0,true);

    CrouzeixRaviartFECollection fec_cr;
    FiniteElementSpace fes_cr(&m, &fec_cr);

    // Boundary conditions: interpolate Dirichlet on boundary edges midpoints
    class DirichletCoeff : public Coefficient {
        const Laplace* L;
    public: explicit DirichletCoeff(const Laplace* l): L(l) {}
        virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) {
            Vector x; T.Transform(ip, x); return Laplace::evalBC(L, x(0), x(1));
        }
    } dbc(this);

    Array<int> ess_tdof_list; // For CR, essential dofs live on boundary faces
    Array<int> ess_bdr;
    if (m.bdr_attributes.Size() == 0) { m.GenerateBoundaryElements(); }
    int nb = m.bdr_attributes.Size() ? m.bdr_attributes.Max() : 1;
    ess_bdr.SetSize(nb); ess_bdr = 1;
    fes_cr.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    BilinearForm a(&fes_cr);
    a.AddDomainIntegrator(new DiffusionIntegrator);
    a.Assemble();

    GridFunction x(&fes_cr); x = 0.0;
    x.ProjectBdrCoefficient(dbc, ess_bdr);

    LinearForm b(&fes_cr); b.Assemble();
    SparseMatrix A; Vector X, B;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
    GSSmoother M(A);
    PCG(A, M, B, X, 1, 2000, 1e-12, 0.0);
    a.RecoverFEMSolution(X, b, x);

    // Transfer CR solution to P1 at vertices: average incident edge dofs via element-local mapping
    u.assign(static_cast<std::size_t>(nv), 0.0);
    std::vector<int> deg(static_cast<std::size_t>(nv), 0);
    for (int el = 0; el < m.GetNE(); ++el) {
        mfem::Array<int> v; m.GetElementVertices(el, v); // v.Size()==3
        mfem::Array<int> edofs; fes_cr.GetElementVDofs(el, edofs); // 3 edge dofs
        // Edge 0: (v0,v1), Edge 1: (v1,v2), Edge 2: (v2,v0)
        int tri_edges[3][2] = { {v[0], v[1]}, {v[1], v[2]}, {v[2], v[0]} };
        for (int e = 0; e < 3 && e < edofs.Size(); ++e) {
            int dof = edofs[e]; if (dof < 0) dof = -1 - dof;
            double val = (dof >= 0 && dof < x.Size()) ? x(dof) : 0.0;
            int v0 = tri_edges[e][0], v1 = tri_edges[e][1];
            u[static_cast<std::size_t>(v0)] += val; deg[static_cast<std::size_t>(v0)]++;
            u[static_cast<std::size_t>(v1)] += val; deg[static_cast<std::size_t>(v1)]++;
        }
    }
    for (int i = 0; i < nv; ++i) { if (deg[static_cast<std::size_t>(i)]>0) u[static_cast<std::size_t>(i)] /= deg[static_cast<std::size_t>(i)]; }
    return true;
}

bool Laplace::writeVTK(const std::string& path) const {
    std::ofstream ofs(path);
    if (!ofs) return false;
    ofs << "# vtk DataFile Version 2.0\n";
    ofs << "xfstream_laplace\n";
    ofs << "ASCII\n";
    ofs << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << nv << " float\n";
    for (int i = 0; i < nv; ++i) {
        ofs << static_cast<float>(verts[static_cast<std::size_t>(i)][0]) << ' '
            << static_cast<float>(verts[static_cast<std::size_t>(i)][1]) << ' ' << 0.0f << "\n";
    }
    // cells: triangles
    const int cellsCount = nelems;
    ofs << "CELLS " << cellsCount << ' ' << cellsCount * 4 << "\n"; // 3 verts + leading count per cell
    for (int e = 0; e < nelems; ++e) {
        auto tri = elems[static_cast<std::size_t>(e)];
        ofs << 3 << ' ' << tri[0] << ' ' << tri[1] << ' ' << tri[2] << "\n";
    }
    ofs << "CELL_TYPES " << nelems << "\n";
    for (int e = 0; e < nelems; ++e) ofs << 5 << "\n"; // VTK_TRIANGLE = 5
    // point data
    ofs << "POINT_DATA " << nv << "\n";
    ofs << "SCALARS u float 1\n";
    ofs << "LOOKUP_TABLE default\n";
    for (int i = 0; i < nv; ++i) ofs << static_cast<float>(u[static_cast<std::size_t>(i)]) << "\n";
    return true;
}
