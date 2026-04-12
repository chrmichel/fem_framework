#pragma once

#include "fem/core/mesh.hpp"
#include "fem/discretization/element/finite_element.hpp"
#include "fem/discretization/quadrature/quadrature_rule.hpp"
#include "fem/linalg/vector.hpp"

#include <cmath>
#include <functional>
#include <vector>

namespace fem::tests {

inline core::Mesh<1> make_uniform_mesh(double a, double b, std::size_t n_elements) {
    std::vector<core::Mesh<1>::Point> nodes;
    nodes.reserve(n_elements + 1);
    const double h = (b - a) / static_cast<double>(n_elements);
    for (std::size_t i = 0; i <= n_elements; ++i)
        nodes.push_back({a + h * static_cast<double>(i)});
    return core::Mesh<1>(std::move(nodes));
}

template<int Dim>
inline double l2_error(
    const fem::core::Mesh<Dim>& mesh,
    const fem::linalg::Vector& u,
    const fem::discretization::element::FiniteElement& fe,
    const fem::discretization::quadrature::QuadratureRule& quad,
    const std::function<double(typename fem::core::Mesh<Dim>::Point)>& u_exact)
{
    double accum = 0.0;
    const std::size_t ndofs = fe.num_dofs();

    for (std::size_t e = 0; e < mesh.n_elements(); ++e) {

        std::vector<std::size_t> g(ndofs);
        for (std::size_t i = 0; i < ndofs; ++i)
            g[i] = mesh.element_node(e, i);

        if constexpr (Dim == 1) {
            const double x0  = mesh.node(g[0])[0];
            const double x1  = mesh.node(g[ndofs - 1])[0];
            const double J   = 0.5 * (x1 - x0);
            const double mid = 0.5 * (x0 + x1);

            for (std::size_t q = 0; q < quad.size(); ++q) {
                const double xi = quad.point(q);
                const double w  = quad.weight(q);
                const double xq = mid + J * xi;

                double uh = 0.0;
                for (std::size_t i = 0; i < ndofs; ++i)
                    uh += u(g[i]) * fe.shape(i, xi);

                const double diff = uh - u_exact({xq});
                accum += diff * diff * (w * J);
            }

        } else if constexpr (Dim == 2) {
            const auto p0 = mesh.node(g[0]);
            const auto p1 = mesh.node(g[1]);
            const auto p2 = mesh.node(g[2]);

            const double J00  = p1[0] - p0[0];
            const double J01  = p2[0] - p0[0];
            const double J10  = p1[1] - p0[1];
            const double J11  = p2[1] - p0[1];
            const double detJ = J00 * J11 - J01 * J10;

            for (std::size_t q = 0; q < quad.size(); ++q) {
                const auto   pt  = quad.point2d(q);
                const double xi  = pt[0];
                const double eta = pt[1];
                const double w   = quad.weight(q);

                double xq = 0.0, yq = 0.0;
                for (std::size_t i = 0; i < ndofs; ++i) {
                    const double phi_i = fe.shape(i, xi, eta);
                    xq += mesh.node(g[i])[0] * phi_i;
                    yq += mesh.node(g[i])[1] * phi_i;
                }

                double uh = 0.0;
                for (std::size_t i = 0; i < ndofs; ++i)
                    uh += u(g[i]) * fe.shape(i, xi, eta);

                const double diff = uh - u_exact({xq, yq});
                accum += diff * diff * (w * detJ);
            }
        }
    }

    return std::sqrt(accum);
}
} // namespace fem::tests