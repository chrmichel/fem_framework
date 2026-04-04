#pragma once
#include "boundary_condition.hpp"
#include <functional>
#include <optional>

namespace fem::boundary {

class NeumannBC final : public BoundaryCondition {
public:
    using Function = std::function<double(double)>;

    NeumannBC(std::optional<Function> g_left,
              std::optional<Function> g_right);

    void apply(linalg::SparseMatrix& A,
               linalg::Vector& b,
               const core::Mesh<1>& mesh) const override;

    std::string name() const override { return "NeumannBC"; }

private:
    std::optional<Function> g_left_;
    std::optional<Function> g_right_;
};

} // namespace fem::boundary