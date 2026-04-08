#pragma once

#include "fem/core/mesh.hpp"
#include "fem/core/mesh_impl.hpp"
#include "fem/problems/problem.hpp"
#include "fem/discretization/element/finite_element.hpp"
#include "fem/discretization/quadrature/quadrature_rule.hpp"
#include "fem/linalg/sparse_matrix.hpp"
#include "fem/linalg/vector.hpp"

#include <stdexcept>
#include <sstream>
#include <vector>

namespace fem::assembly {

template<int Dim>
class Assembler {
public:
    Assembler(const core::Mesh<Dim>& mesh,
              const problems::Problem& problem,
              const discretization::element::FiniteElement& fe,
              const discretization::quadrature::QuadratureRule& quad)
        : m_mesh(mesh),
          m_problem(problem),
          m_fe(fe),
          m_quad(quad)
    {}

    void assemble(linalg::SparseMatrix& A,
              linalg::Vector& b) const
    {
        const std::size_t ndofs_local = m_fe.num_dofs();

        for (std::size_t e = 0; e < m_mesh.n_elements(); ++e) {

            std::vector<std::size_t> g(ndofs_local);
            for (std::size_t i = 0; i < ndofs_local; ++i)
                g[i] = m_mesh.element_node(e, i);

            std::vector<double> fe(ndofs_local, 0.0);
            std::vector<double> phi(ndofs_local, 0.0);
            std::vector<double> dphidx(ndofs_local, 0.0);
            std::vector<std::vector<double>> Ke(
                ndofs_local, std::vector<double>(ndofs_local, 0.0));

            if constexpr (Dim == 1) {

                const double x0 = m_mesh.node(g[0])[0];
                const double x1 = m_mesh.node(g[ndofs_local - 1])[0];
                const double J  = 0.5 * (x1 - x0);

                if (J <= 0.0) {
                    std::ostringstream oss;
                    oss << "Assembler: non-positive Jacobian on element " << e;
                    throw std::runtime_error(oss.str());
                }

                for (std::size_t q = 0; q < m_quad.size(); ++q) {
                    const double xi = m_quad.point(q);
                    const double w  = m_quad.weight(q);

                    double xq = 0.0;
                    for (std::size_t i = 0; i < ndofs_local; ++i)
                        xq += m_mesh.node(g[i])[0] * m_fe.shape(i, xi);

                    const double a  = m_problem.diffusion(xq);
                    const double fq = m_problem.rhs(xq);
                    const double c  = m_problem.reaction(xq);

                    for (std::size_t i = 0; i < ndofs_local; ++i) {
                        phi[i]    = m_fe.shape(i, xi);
                        dphidx[i] = m_fe.dshape_dxi(i, xi) / J;
                    }

                    const double weight_phys = w * J;

                    for (std::size_t i = 0; i < ndofs_local; ++i) {
                        fe[i] += fq * phi[i] * weight_phys;
                        for (std::size_t j = 0; j < ndofs_local; ++j) {
                            Ke[i][j] += a * dphidx[i] * dphidx[j] * weight_phys;
                            Ke[i][j] += c * phi[i]    * phi[j]    * weight_phys;
                        }
                    }
                }

            } else if constexpr (Dim == 2) {

                std::vector<double> dphidy(ndofs_local, 0.0);

                const auto p0 = m_mesh.node(g[0]);
                const auto p1 = m_mesh.node(g[1]);
                const auto p2 = m_mesh.node(g[2]);

                const double J00 = p1[0] - p0[0];
                const double J01 = p2[0] - p0[0];
                const double J10 = p1[1] - p0[1];
                const double J11 = p2[1] - p0[1];

                const double detJ = J00 * J11 - J01 * J10;

                if (detJ <= 0.0) {
                    std::ostringstream oss;
                    oss << "Assembler: non-positive Jacobian determinant on element " << e;
                    throw std::runtime_error(oss.str());
                }

                const double invJ00 =  J11 / detJ;
                const double invJ01 = -J01 / detJ;
                const double invJ10 = -J10 / detJ;
                const double invJ11 =  J00 / detJ;

                for (std::size_t q = 0; q < m_quad.size(); ++q) {
                    const auto   pt  = m_quad.point2d(q);
                    const double xi  = pt[0];
                    const double eta = pt[1];
                    const double w   = m_quad.weight(q);

                    double xq = 0.0, yq = 0.0;
                    for (std::size_t i = 0; i < ndofs_local; ++i) {
                        const double phi_i = m_fe.shape(i, xi, eta);
                        xq += m_mesh.node(g[i])[0] * phi_i;
                        yq += m_mesh.node(g[i])[1] * phi_i;
                    }

                    const double a  = m_problem.diffusion(xq, yq);
                    const double fq = m_problem.rhs(xq, yq);
                    const double c  = m_problem.reaction(xq, yq);

                    for (std::size_t i = 0; i < ndofs_local; ++i) {
                        phi[i] = m_fe.shape(i, xi, eta);

                        const double dphi_dxi  = m_fe.dshape_dxi(i, xi, eta);
                        const double dphi_deta = m_fe.dshape_deta(i, xi, eta);

                        dphidx[i] = invJ00 * dphi_dxi + invJ10 * dphi_deta;
                        dphidy[i] = invJ01 * dphi_dxi + invJ11 * dphi_deta;
                    }

                    const double weight_phys = w * detJ;

                    for (std::size_t i = 0; i < ndofs_local; ++i) {
                        fe[i] += fq * phi[i] * weight_phys;
                        for (std::size_t j = 0; j < ndofs_local; ++j) {
                            Ke[i][j] += a * (dphidx[i] * dphidx[j] +
                                            dphidy[i] * dphidy[j]) * weight_phys;
                            Ke[i][j] += c * phi[i] * phi[j] * weight_phys;
                        }
                    }
                }
            }

            // Scatter — gemeinsam für 1D und 2D
            for (std::size_t i = 0; i < ndofs_local; ++i) {
                b(g[i]) += fe[i];
                for (std::size_t j = 0; j < ndofs_local; ++j)
                    A(g[i], g[j]) += Ke[i][j];
            }
        }
    }

private:
    const core::Mesh<Dim>& m_mesh;
    const problems::Problem& m_problem;
    const discretization::element::FiniteElement& m_fe;
    const discretization::quadrature::QuadratureRule& m_quad;
};

extern template class Assembler<1>;
extern template class Assembler<2>;
extern template class Assembler<3>;

} // namespace fem::assembly