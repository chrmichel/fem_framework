#pragma once

#include "fem/core/mesh.hpp"
#include "fem/problems/problem.hpp"
#include "fem/discretization/element/finite_element.hpp"
#include "fem/discretization/quadrature/quadrature_rule.hpp"
#include "fem/linalg/matrix.hpp"
#include "fem/linalg/vector.hpp"

namespace fem::assembly {

class Assembler {
/**
 * Generic FEM Assembler (Phase 2).
 *
 * Works with:
 *   - Mesh (geometry)
 *   - Problem (physics)
 *   - FiniteElement (reference basis)
 *   - QuadratureRule (numerical integration)
 *
 * No hardcoded P1 logic allowed anymore.
 */
public:
    Assembler(const core::Mesh& mesh,
              const problems::Problem& problem,
              const discretization::element::FiniteElement& fe,
              const discretization::quadrature::QuadratureRule& quad);

    /// Assemble global stiffness matrix and RHS vector
    void assemble(linalg::Matrix& A,
                  linalg::Vector& b) const;

private:
    const core::Mesh& m_mesh;
    const problems::Problem& m_problem;
    const discretization::element::FiniteElement& m_fe;
    const discretization::quadrature::QuadratureRule& m_quad;
};

} // namespace fem::assembly