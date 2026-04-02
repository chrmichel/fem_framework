#include "fem/boundary/neumann_bc.hpp"
#include "fem/core/mesh.hpp"
#include "fem/linalg/vector.hpp"
#include "fem/linalg/matrix.hpp"

#include <stdexcept>
#include <utility>

namespace fem::boundary {

NeumannBC::NeumannBC(std::optional<Function> g_left,
                     std::optional<Function> g_right)
    : g_left_(std::move(g_left)), g_right_(std::move(g_right))
{
    if (!g_left_ && !g_right_) {
        throw std::invalid_argument("NeumannBC: at least one boundary function must be set.");
    }
}

void NeumannBC::apply(linalg::Matrix& /*A*/,
                      linalg::Vector& b,
                      const core::Mesh& mesh) const
{
    if (g_left_) {
        const std::size_t iL = mesh.left_boundary_node();
        b(iL) -= (*g_left_)(mesh.node(iL));
    }

    if (g_right_) {
        const std::size_t iR = mesh.right_boundary_node();
        b(iR) += (*g_right_)(mesh.node(iR));
    }
}

} // namespace fem::boundary