#include "fem/assembly/assembler.hpp"

#include <stdexcept>
#include <sstream>

namespace fem::assembly {

Assembler::Assembler(const core::Mesh& mesh,
                     const problems::Problem& problem,
                     const discretization::element::FiniteElement& fe,
                     const discretization::quadrature::QuadratureRule& quad)
    : m_mesh(mesh),
      m_problem(problem),
      m_fe(fe),
      m_quad(quad)
{
}

void Assembler::assemble(linalg::Matrix& A,
                         linalg::Vector& b) const
{
    const std::size_t ndofs_local = m_fe.num_dofs();

    // Phase-2 Constraint:
    // Mesh connectivity is still P1-based (2 nodes per element).
    if (ndofs_local != 2) {
        std::ostringstream oss;
        oss << "Assembler: Mesh connectivity is P1-only but FiniteElement has "
            << ndofs_local << " local DoFs.";
        throw std::invalid_argument(oss.str());
    }

    for (std::size_t e = 0; e < m_mesh.n_elements(); ++e) {

        const std::size_t g0 = m_mesh.element_node(e, 0);
        const std::size_t g1 = m_mesh.element_node(e, 1);

        const double x0 = m_mesh.node(g0);
        const double x1 = m_mesh.node(g1);

        const double J = 0.5 * (x1 - x0);
        if (J <= 0.0) {
            std::ostringstream oss;
            oss << "Assembler: non-positive Jacobian on element " << e;
            throw std::runtime_error(oss.str());
        }

        double Ke[2][2] = {{0.0, 0.0},
                           {0.0, 0.0}};
        double fe[2]    = {0.0, 0.0};

        // Quadrature loop
        for (std::size_t q = 0; q < m_quad.size(); ++q) {

            const double xi = m_quad.point(q);
            const double w  = m_quad.weight(q);

            const double xq = 0.5 * (x0 + x1) + J * xi;

            const double a  = m_problem.diffusion(xq);
            const double fq = m_problem.rhs(xq);
            const double c  = m_problem.reaction(xq);   // NEW

            double phi[2];
            double dphidx[2];

            for (std::size_t i = 0; i < 2; ++i) {
                phi[i] = m_fe.shape(i, xi);
                const double dphi_dxi = m_fe.dshape_dxi(i, xi);
                dphidx[i] = dphi_dxi / J;
            }

            const double weight_phys = w * J;

            for (std::size_t i = 0; i < 2; ++i) {
                fe[i] += fq * phi[i] * weight_phys;

                for (std::size_t j = 0; j < 2; ++j) {
                    Ke[i][j] += a * dphidx[i] * dphidx[j] * weight_phys;
                    // reaction term: c * u v   (NEW)
                    Ke[i][j] += c * phi[i] * phi[j] * weight_phys;
                }
            }
        }

        const std::size_t g[2] = {g0, g1};

        for (std::size_t i = 0; i < 2; ++i) {
            b(g[i]) += fe[i];

            for (std::size_t j = 0; j < 2; ++j) {
                A(g[i], g[j]) += Ke[i][j];
            }
        }
    }
}

} // namespace fem::assembly