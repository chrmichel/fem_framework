#include "fem/analysis/refine.hpp"

#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace fem::analysis {

core::Mesh1D refine(const core::Mesh1D& mesh,
                  const ErrorEstimate& estimate,
                  double theta)
{
    if (estimate.indicators.size() != mesh.n_elements()) {
        throw std::invalid_argument(
            "refine: indicators size does not match mesh.n_elements()");
    }

    if (theta <= 0.0 || theta > 1.0) {
        throw std::invalid_argument(
            "refine: theta must be in (0, 1]");
    }

    const double max_eta = *std::max_element(
        estimate.indicators.begin(), estimate.indicators.end());

    if (max_eta == 0.0) {
        return mesh;
    }

    const double threshold = theta * max_eta;

    std::vector<core::Mesh1D::Point> new_nodes;
    std::vector<std::vector<std::size_t>> new_conn;

    new_nodes.reserve(mesh.n_nodes());
    new_conn.reserve(mesh.n_elements());

    for (std::size_t i = 0; i < mesh.n_nodes(); ++i)
        new_nodes.push_back(mesh.node(i));

    for (std::size_t e = 0; e < mesh.n_elements(); ++e) {
        const std::size_t g0   = mesh.element_node(e, 0);
        const std::size_t g1   = mesh.element_node(e, mesh.nodes_per_element(e) - 1);

        if (estimate.indicators[e] > threshold) {
            const double xmid    = 0.5 * (mesh.node(g0)[0] + mesh.node(g1)[0]);
            const std::size_t gm = new_nodes.size();
            new_nodes.push_back({xmid});

            new_conn.push_back({g0, gm});
            new_conn.push_back({gm, g1});
        } else {
            new_conn.push_back({g0, g1});
        }
    }

    for (std::size_t i = 0; i < new_nodes.size(); ++i)
        std::cerr << "node " << i << " = " << new_nodes[i][0] << "\n";
    for (std::size_t e = 0; e < new_conn.size(); ++e) {
        std::cerr << "elem " << e << " -> ";
        for (auto idx : new_conn[e]) std::cerr << idx << " ";
        std::cerr << "\n";
    }
    core::Mesh1D refined(std::move(new_nodes), std::move(new_conn));
    refined.set_boundary_nodes(mesh.left_boundary_node(), mesh.right_boundary_node());
    return refined;
}

} // namespace fem::analysis