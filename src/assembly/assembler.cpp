#include "fem/assembly/assembler.hpp"

#include "fem/core/mesh.hpp"
#include "fem/boundary/boundary_condition.hpp"
#include <array>
#include <stdexcept>

namespace fem::assembly {

double Assembler::phi(int a, double xi)
{
    if (a == 0) return 0.5 * (1.0 - xi);
    if (a == 1) return 0.5 * (1.0 + xi);
    throw std::invalid_argument("Assembler::phi: invalid local dof index");
}

double Assembler::dphi_dxi(int a)
{
    if (a == 0) return -0.5;
    if (a == 1) return  0.5;
    throw std::invalid_argument("Assembler::dphi_dxi: invalid local dof index");
}

fem::linalg::Matrix Assembler::assemble_stiffness(const fem::problems::Problem& problem) const
{
    const fem::Mesh& mesh = problem.mesh();

    if (mesh.dimension() != 1) {
        throw std::invalid_argument("Assembler: Phase 1 supports only 1D meshes.");
    }

    const std::size_t N = mesh.n_nodes();
    fem::linalg::Matrix A(N, N, 0.0);

    // 2-Punkt Gauss auf [-1,1]
    constexpr std::array<double,2> xi = { -0.5773502691896257, 0.5773502691896257 };
    constexpr std::array<double,2> w  = {  1.0,                 1.0 };

    for (std::size_t e = 0; e < mesh.n_elements(); ++e) {
        const auto conn = mesh.element(e); // {i0,i1}
        const std::size_t i0 = conn[0];
        const std::size_t i1 = conn[1];

        const double x0 = mesh.node(i0);
        const double x1 = mesh.node(i1);
        const double J  = (x1 - x0) / 2.0; // dx/dxi

        if (J <= 0.0) {
            throw std::runtime_error("Assembler: non-positive Jacobian.");
        }

        double Ke[2][2] = { {0.0, 0.0}, {0.0, 0.0} };

        for (int q = 0; q < 2; ++q) {
            const double x = 0.5 * (x0 + x1) + J * xi[q];
            const double a = problem.diffusion(x);

            const double dN_dx[2] = {
                dphi_dxi(0) / J,
                dphi_dxi(1) / J
            };

            for (int Aidx = 0; Aidx < 2; ++Aidx) {
                for (int Bidx = 0; Bidx < 2; ++Bidx) {
                    Ke[Aidx][Bidx] += w[q] * J * a * dN_dx[Aidx] * dN_dx[Bidx];
                }
            }
        }

        const std::size_t gdof[2] = { i0, i1 };
        for (int a = 0; a < 2; ++a) {
            for (int b = 0; b < 2; ++b) {
                A(gdof[a], gdof[b]) += Ke[a][b];
            }
        }
    }

    return A;
}

fem::linalg::Vector Assembler::assemble_load(const fem::problems::Problem& problem) const
{
    const fem::Mesh& mesh = problem.mesh();

    if (mesh.dimension() != 1) {
        throw std::invalid_argument("Assembler: Phase 1 supports only 1D meshes.");
    }

    const std::size_t N = mesh.n_nodes();
    fem::linalg::Vector b(N, 0.0);

    // 2-Punkt Gauss auf [-1,1]
    constexpr std::array<double,2> xi = { -0.5773502691896257, 0.5773502691896257 };
    constexpr std::array<double,2> w  = {  1.0,                 1.0 };

    for (std::size_t e = 0; e < mesh.n_elements(); ++e) {
        const auto conn = mesh.element(e);
        const std::size_t i0 = conn[0];
        const std::size_t i1 = conn[1];

        const double x0 = mesh.node(i0);
        const double x1 = mesh.node(i1);
        const double J  = (x1 - x0) / 2.0;

        if (J <= 0.0) {
            throw std::runtime_error("Assembler: non-positive Jacobian.");
        }

        double Fe[2] = {0.0, 0.0};

        for (int q = 0; q < 2; ++q) {
            const double x = 0.5 * (x0 + x1) + J * xi[q];
            const double f = problem.rhs(x);

            const double Nref[2] = { phi(0, xi[q]), phi(1, xi[q]) };

            for (int a = 0; a < 2; ++a) {
                Fe[a] += w[q] * J * f * Nref[a];
            }
        }

        b(i0) += Fe[0];
        b(i1) += Fe[1];
    }

    return b;
}

void Assembler::assemble_system(const fem::problems::Problem& problem,
                               fem::linalg::Matrix& A,
                               fem::linalg::Vector& b) const
{
    A = assemble_stiffness(problem);
    b = assemble_load(problem);

    // BCs werden zentral über das Problem geliefert
    for (const auto& bc : problem.boundary_conditions()) {
        bc->apply(A, b, problem.mesh());
    }
}

} // namespace fem::assembly