#pragma once
#include "boundary_condition.hpp"

#include "fem/core/mesh.hpp"
#include "fem/linalg/sparse_matrix.hpp"
#include "fem/linalg/vector.hpp"
#include "fem/linalg/checks.hpp"
#include <functional>
#include <optional>

namespace fem::boundary {

template<int Dim>
class NeumannBC final : public BoundaryCondition<Dim> {
public:
    using Point    = typename core::Mesh<Dim>::Point;
    using Function = std::function<double(Point)>;

    NeumannBC(std::optional<Function> g_left,
              std::optional<Function> g_right)
        : g_left_(std::move(g_left)), g_right_(std::move(g_right))
    {
        if (!g_left_ && !g_right_) {
            throw std::invalid_argument("NeumannBC: at least one boundary function must be set.");
        }
    }

    void apply(linalg::SparseMatrix& A,
               linalg::Vector& b,
               const core::Mesh<Dim>& mesh) const override
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

    std::string name() const override { return "NeumannBC"; }

private:
    std::optional<Function> g_left_;
    std::optional<Function> g_right_;
};

} // namespace fem::boundary