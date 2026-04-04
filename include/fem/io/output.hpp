#pragma once

#include "fem/core/mesh.hpp"
#include "fem/problems/problem.hpp"
#include "fem/discretization/element/finite_element.hpp"
#include "fem/discretization/quadrature/quadrature_rule.hpp"
#include "fem/boundary/boundary_condition.hpp"
#include "fem/linalg/vector.hpp"

#include <filesystem>
#include <string>

namespace fem::io {

std::filesystem::path write_results(
    const core::Mesh1D& mesh,
    const linalg::Vector& solution,
    const problems::Problem& problem,
    const discretization::element::FiniteElement& fe,
    const discretization::quadrature::QuadratureRule& quad,
    const std::vector<boundary::BoundaryCondition*>& bcs,
    const std::string& base_dir = "results");

} // namespace fem::io