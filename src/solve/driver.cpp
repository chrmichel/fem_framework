#include "fem/solve/driver.hpp"

#include "fem/assembly/assembler.hpp"
#include "fem/linalg/matrix.hpp"
#include "fem/linalg/vector.hpp"
#include "fem/linalg/solver.hpp"   // <- falls deine freie Funktion hier deklariert ist

namespace fem {

linalg::Vector Driver::solve(
    const core::Mesh& mesh,
    const problems::Problem& problem,
    const discretization::element::FiniteElement& fe,
    const discretization::quadrature::QuadratureRule& quad,
    const boundary::BoundaryCondition& bc)
{
    linalg::Matrix A(mesh.n_nodes(), mesh.n_nodes());
    linalg::Vector b(mesh.n_nodes());

    assembly::Assembler assembler(mesh, problem, fe, quad);
    assembler.assemble(A, b);

    bc.apply(A, b, mesh);

    // Free function solver
    return linalg::solve(A, b);
}

} // namespace fem