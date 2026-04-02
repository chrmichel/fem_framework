#pragma once

#include "fem/core/mesh.hpp"

namespace fem::core {
    Mesh create_p1_mesh_1d(double a, double b, std::size_t num_elements);
    Mesh create_p2_mesh_1d(double a, double b, std::size_t num_elements);

} // namespace fem::core