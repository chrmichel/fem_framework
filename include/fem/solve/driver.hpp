#pragma once

#include <vector>

#include "fem/core/mesh.hpp"
#include "fem/core/mesh_impl.hpp"
#include "fem/problems/problem.hpp"
#include "fem/discretization/element/finite_element.hpp"
#include "fem/discretization/quadrature/quadrature_rule.hpp"
#include "fem/boundary/boundary_condition.hpp"
#include "fem/linalg/vector.hpp"
#include "fem/linalg/sparse_matrix.hpp"
#include "fem/linalg/solver.hpp"
#include "fem/assembly/assembler.hpp"
#include "fem/assembly/sparsity.hpp"

namespace fem {

template<int Dim>
class Driver {
public:
    static linalg::Vector solve(
        const core::Mesh<Dim>& mesh,
        const problems::Problem& problem,
        const discretization::element::FiniteElement& fe,
        const discretization::quadrature::QuadratureRule& quad,
        const std::vector<const boundary::BoundaryCondition<Dim>*>& bcs)
    {
        linalg::SparseMatrix A = assembly::build_sparse_matrix(mesh, fe);
        linalg::Vector b(mesh.n_nodes());

        assembly::Assembler<Dim> assembler(mesh, problem, fe, quad);
        assembler.assemble(A, b);

        for (const auto* bc : bcs)
            bc->apply(A, b, mesh);

        return linalg::solve(A, b);
    }
};

} // namespace fem