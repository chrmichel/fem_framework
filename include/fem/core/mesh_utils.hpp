#pragma once

#include "fem/core/mesh.hpp"

namespace fem::core {
    Mesh1D create_p1_mesh_1d(double a, double b, std::size_t num_elements);
    Mesh1D create_p2_mesh_1d(double a, double b, std::size_t num_elements);

    Mesh<2> create_triangle_mesh_2d(
    double x0, double x1,
    double y0, double y1,
    std::size_t nx, std::size_t ny);

} // namespace fem::core