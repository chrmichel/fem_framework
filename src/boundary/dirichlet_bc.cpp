#include "fem/boundary/dirichlet_bc.hpp"

#include "fem/core/mesh.hpp"
#include "fem/linalg/matrix.hpp"
#include "fem/linalg/vector.hpp"
#include "fem/linalg/checks.hpp"

#include <stdexcept>
#include <utility>

namespace fem::boundary {

DirichletBC::DirichletBC(Function g_left, Function g_right)
    : g_left_(std::move(g_left)), g_right_(std::move(g_right))
{
    if (!g_left_ || !g_right_) {
        throw std::invalid_argument("DirichletBC: boundary functions must be valid.");
    }
}

void DirichletBC::apply(linalg::Matrix& A,
                        linalg::Vector& b,
                        const core::Mesh& mesh) const
{
    // Checks live in linalg layer (as requested)
    linalg::require_square(A);
    linalg::require_same_size(A, b);

    const std::size_t iL = mesh.left_boundary_node();
    const std::size_t iR = mesh.right_boundary_node();

    // left boundary value
    {
        const double xL = mesh.node(iL);
        const double vL = g_left_(xL);
        enforce_node(iL, vL, A, b);
    }

    // right boundary value (avoid double apply for degenerate meshes)
    if (iR != iL) {
        const double xR = mesh.node(iR);
        const double vR = g_right_(xR);
        enforce_node(iR, vR, A, b);
    }
}

void DirichletBC::enforce_node(std::size_t node,
                               double value,
                               linalg::Matrix& A,
                               linalg::Vector& b)
{
    // Inhomogeneous Dirichlet elimination with RHS correction:
    // For all j != node: b_j -= A(j,node) * value, then A(j,node) = 0
    const std::size_t n = A.rows();

    for (std::size_t j = 0; j < n; ++j) {
        if (j == node) continue;
        b(j) -= A(j, node) * value;
        A(j, node) = 0.0;
    }

    // Zero out row 'node'
    A.set_row_zero(node);

    // Set diagonal and rhs
    A(node, node) = 1.0;
    b(node) = value;
}

} // namespace fem::boundary