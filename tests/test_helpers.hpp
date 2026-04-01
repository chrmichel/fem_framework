#pragma once

#include "fem/core/mesh.hpp"
#include "fem/discretization/element/finite_element.hpp"
#include "fem/discretization/quadrature/quadrature_rule.hpp"
#include "fem/linalg/vector.hpp"

#include <cmath>
#include <functional>
#include <vector>

namespace fem::tests {

inline core::Mesh make_uniform_mesh(double a, double b, std::size_t n_elements) {
    std::vector<double> nodes;
    nodes.reserve(n_elements + 1);
    const double h = (b - a) / static_cast<double>(n_elements);
    for (std::size_t i = 0; i <= n_elements; ++i) {
        nodes.push_back(a + h * static_cast<double>(i));
    }
    return core::Mesh(std::move(nodes));
}

} // namespace fem::tests