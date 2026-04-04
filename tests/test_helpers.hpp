#pragma once

#include "fem/core/mesh.hpp"
#include "fem/discretization/element/finite_element.hpp"
#include "fem/discretization/quadrature/quadrature_rule.hpp"
#include "fem/linalg/vector.hpp"

#include <cmath>
#include <functional>
#include <vector>

namespace fem::tests {

inline core::Mesh1D make_uniform_mesh(double a, double b, std::size_t n_elements) {
    std::vector<core::Mesh1D::Point> nodes;
    nodes.reserve(n_elements + 1);
    const double h = (b - a) / static_cast<double>(n_elements);
    for (std::size_t i = 0; i <= n_elements; ++i) {
        nodes.push_back({a + h * static_cast<double>(i)});
    }
    return core::Mesh1D(std::move(nodes));
}

inline double l2_error(const fem::core::Mesh1D& mesh,
                       const fem::linalg::Vector& u,
                       const fem::discretization::element::FiniteElement& fe,
                       const fem::discretization::quadrature::QuadratureRule& quad,
                       const std::function<double(double)>& u_exact)
{
    double accum = 0.0;
    const std::size_t ndofs = fe.num_dofs();

    for (std::size_t e = 0; e < mesh.n_elements(); ++e) {

        const std::size_t g0 = mesh.element_node(e, 0);
        const std::size_t glast = mesh.element_node(e, ndofs - 1);

        const double x0 = mesh.node(g0)[0];
        const double x1 = mesh.node(glast)[0];
        const double J   = 0.5 * (x1 - x0);
        const double mid = 0.5 * (x0 + x1);

        for (std::size_t q = 0; q < quad.size(); ++q) {
            const double xi = quad.point(q);
            const double w  = quad.weight(q);

            const double xq = mid + J * xi;

            double uh = 0.0;
            for (std::size_t i = 0; i < ndofs; ++i)
                uh += u(mesh.element_node(e, i)) * fe.shape(i, xi);

            const double diff = uh - u_exact(xq);
            accum += diff * diff * (w * J);
        }
    }

    return std::sqrt(accum);
}

} // namespace fem::tests