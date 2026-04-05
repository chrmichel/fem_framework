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

            const double x0 = m_mesh.node(g[0])[0];
            const double x1 = m_mesh.node(g[ndofs_local - 1])[0];

            const double J = 0.5 * (x1 - x0);
            if (J <= 0.0) {
                std::ostringstream oss;
                oss << "Assembler: non-positive Jacobian on element " << e;
                throw std::runtime_error(oss.str());
            }

            std::vector<double> fe(ndofs_local, 0.0);
            std::vector<double> phi(ndofs_local, 0.0);
            std::vector<double> dphidx(ndofs_local, 0.0);
            std::vector<std::vector<double>> Ke(
                ndofs_local, std::vector<double>(ndofs_local, 0.0));

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