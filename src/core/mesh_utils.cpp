#include "fem/core/mesh_utils.hpp"

namespace fem::core {
Mesh1D create_p1_mesh_1d(double a, double b, std::size_t num_elements) {
    std::vector<Mesh1D::Point> nodes;
    const double h = (b - a) / num_elements;
    for (std::size_t i = 0; i <= num_elements; ++i) {
        nodes.push_back({a + i * h});
    }
    return Mesh1D(nodes);
}

Mesh1D create_p2_mesh_1d(double a, double b, std::size_t num_elements) {
    std::vector<Mesh1D::Point> nodes;
    const double h = (b - a) / num_elements / 2.0;
    for (std::size_t i = 0; i <= 2 * num_elements; ++i) {
        nodes.push_back({a + i * h});
    }
    std::vector<std::vector<std::size_t>> connectivity;
    for (std::size_t e = 0; e < num_elements; ++e) {
        connectivity.push_back({2*e, 2*e+1, 2*e+2});
    }
    Mesh1D m(std::move(nodes), std::move(connectivity));
    m.set_boundary_nodes(0, 2*num_elements); // first and last node are boundaries
    return m;
}

Mesh<2> create_triangle_mesh_2d(
    double x0, double x1,
    double y0, double y1,
    std::size_t nx, std::size_t ny)
{
    using Point = Mesh<2>::Point;

    // Knoten aufbauen
    std::vector<Point> nodes;
    nodes.reserve((nx + 1) * (ny + 1));

    const double hx = (x1 - x0) / nx;
    const double hy = (y1 - y0) / ny;

    for (std::size_t j = 0; j <= ny; ++j)
        for (std::size_t i = 0; i <= nx; ++i)
            nodes.push_back({x0 + i * hx, y0 + j * hy});

    // Konnektivität aufbauen — 2 Dreiecke pro Rechteck
    std::vector<std::vector<std::size_t>> connectivity;
    connectivity.reserve(2 * nx * ny);

    auto idx = [&](std::size_t i, std::size_t j) {
        return j * (nx + 1) + i;
    };

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            // Dreieck 1: untere rechte Hälfte
            connectivity.push_back({idx(i,j), idx(i+1,j), idx(i+1,j+1)});
            // Dreieck 2: obere linke Hälfte
            connectivity.push_back({idx(i,j), idx(i+1,j+1), idx(i,j+1)});
        }
    }

    Mesh<2> mesh(std::move(nodes), std::move(connectivity));
    // Randknoten: alle Knoten auf den vier Seiten des Rechtecks
    std::vector<std::size_t> bnd_ids;
    for (std::size_t j = 0; j <= ny; ++j) {
        for (std::size_t i = 0; i <= nx; ++i) {
            if (i == 0 || i == nx || j == 0 || j == ny)
                bnd_ids.push_back(idx(i, j));
        }
    }
    mesh.set_boundary_node_ids(std::move(bnd_ids));
    return mesh;
}

} // namespace fem::core