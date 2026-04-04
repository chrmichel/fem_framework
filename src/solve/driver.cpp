#include "fem/solve/driver.hpp"

#include "fem/assembly/assembler.hpp"
#include "fem/assembly/sparsity.hpp"
#include "fem/linalg/vector.hpp"
#include "fem/linalg/solver.hpp"   // <- falls deine freie Funktion hier deklariert ist

namespace fem {

linalg::Vector Driver::solve(
    const core::Mesh1D& mesh,
    const problems::Problem& problem,
    const discretization::element::FiniteElement& fe,
    const discretization::quadrature::QuadratureRule& quad,
    const std::vector<boundary::BoundaryCondition*>& bcs)
{
    linalg::SparseMatrix A = assembly::build_sparse_matrix(mesh, fe);
    linalg::Vector b(mesh.n_nodes());

    assembly::Assembler assembler(mesh, problem, fe, quad);
    assembler.assemble(A, b);

    for (const auto* bc : bcs)
        bc->apply(A, b, mesh);

    return linalg::solve(A, b);
}

} // namespace fem