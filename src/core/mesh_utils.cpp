#include "fem/core/mesh_utils.hpp"

namespace fem::core {
Mesh create_p1_mesh_1d(double a, double b, std::size_t num_elements) {
    std::vector<double> nodes;
    const double h = (b - a) / num_elements;
    for (std::size_t i = 0; i <= num_elements; ++i) {
        nodes.push_back(a + i * h);
    }
    return Mesh(nodes);
}

Mesh create_p2_mesh_1d(double a, double b, std::size_t num_elements) {
    std::vector<double> nodes;
    const double h = (b - a) / num_elements / 2.0;
    for (std::size_t i = 0; i <= 2 * num_elements; ++i) {
        double x = a + i * h;
        nodes.push_back(x);
    }
    std::vector<std::vector<std::size_t>> connectivity;
    for (std::size_t e = 0; e < num_elements; ++e) {
        connectivity.push_back({2*e, 2*e+1, 2*e+2});
    }
    Mesh m(std::move(nodes), std::move(connectivity));
    m.set_boundary_nodes(0, 2*num_elements); // first and last node are boundaries
    return m;
}
} // namespace fem::core