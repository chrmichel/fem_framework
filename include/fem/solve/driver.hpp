#pragma once

#include<vector>

#include "fem/core/mesh.hpp"
#include "fem/problems/problem.hpp"
#include "fem/discretization/element/finite_element.hpp"
#include "fem/discretization/quadrature/quadrature_rule.hpp"
#include "fem/boundary/boundary_condition.hpp"
#include "fem/linalg/vector.hpp"

namespace fem {

/**
 * High-level facade for solving FEM systems.
 *
 * Phase 2:
 *   Injection of FiniteElement and QuadratureRule.
 */
class Driver {
public:
    static linalg::Vector solve(
        const core::Mesh& mesh,
        const problems::Problem& problem,
        const discretization::element::FiniteElement& fe,
        const discretization::quadrature::QuadratureRule& quad,
        const std::vector<boundary::BoundaryCondition*>& bcs);
};

} // namespace fem