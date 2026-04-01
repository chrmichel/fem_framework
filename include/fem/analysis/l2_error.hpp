#pragma once
#include "fem/core/mesh.hpp"
#include "fem/discretization/element/finite_element.hpp"
#include "fem/discretization/quadrature/quadrature_rule.hpp"
#include "fem/linalg/vector.hpp"
#include <cmath>
#include <functional>

namespace fem::analysis {
    
inline double l2_error(const fem::core::Mesh& mesh,
                       const fem::linalg::Vector& u,
                       const fem::discretization::element::FiniteElement& fe,
                       const fem::discretization::quadrature::QuadratureRule& quad,
                       const std::function<double(double)>& u_exact)
{
    double accum = 0.0;

    for (std::size_t e = 0; e < mesh.n_elements(); ++e) {
        const std::size_t g0 = mesh.element_node(e, 0);
        const std::size_t g1 = mesh.element_node(e, 1);

        const double x0 = mesh.node(g0);
        const double x1 = mesh.node(g1);

        const double J = 0.5 * (x1 - x0);
        const double mid = 0.5 * (x0 + x1);

        const double u0 = u(g0);
        const double u1 = u(g1);

        for (std::size_t q = 0; q < quad.size(); ++q) {
            const double xi = quad.point(q);
            const double w  = quad.weight(q);

            const double xq = mid + J * xi;

            const double uh = u0 * fe.shape(0, xi) + u1 * fe.shape(1, xi);
            const double diff = uh - u_exact(xq);

            accum += diff * diff * (w * J);
        }
    }

    return std::sqrt(accum);
}
} // namespace fem::analysis
