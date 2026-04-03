#include "fem/analysis/jump_indicator.hpp"

#include <cmath>
#include <numeric>

namespace fem::analysis {

ErrorEstimate jump_indicator(
    const core::Mesh& mesh,
    const linalg::Vector& u,
    const problems::Problem& problem)
{
    const std::size_t n = mesh.n_elements();

    ErrorEstimate est;
    est.indicators.resize(n, 0.0);

    for (std::size_t e = 0; e < n - 1; ++e) {
        const std::size_t g0 = mesh.element_node(e, 0);
        const std::size_t g1 = mesh.element_node(e, mesh.nodes_per_element(e) - 1);
        const std::size_t g2 = mesh.element_node(e + 1, mesh.nodes_per_element(e + 1) - 1);

        const double x0 = mesh.node(g0);
        const double x1 = mesh.node(g1);
        const double x2 = mesh.node(g2);

        const double h_e    = x1 - x0;
        const double h_next = x2 - x1;

        const double du_e    = (u(g1) - u(g0)) / h_e;
        const double du_next = (u(g2) - u(g1)) / h_next;

        const double a_e    = problem.diffusion(0.5 * (x0 + x1));
        const double a_next = problem.diffusion(0.5 * (x1 + x2));

        const double jump      = std::abs(a_e * du_e - a_next * du_next);
        est.indicators[e]      = h_e * jump;
    }

    double global = 0.0;
    for (double eta : est.indicators)
        global += eta * eta;
    est.global_error = std::sqrt(global);

    return est;
}

} // namespace fem::analysis