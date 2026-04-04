#pragma once

#include "boundary_condition.hpp"
#include <functional>
#include <cstddef>
#include <optional>

namespace fem {
class Mesh;

namespace boundary {

class DirichletBC final : public BoundaryCondition {
public:
    using Function = std::function<double(double)>;

    DirichletBC(std::optional<Function> g_left, std::optional<Function> g_right);

    void apply(linalg::SparseMatrix& A,
               linalg::Vector& b,
               const core::Mesh<1>& mesh) const override;
    
    virtual std::string name() const override {
        return "DirichletBC";
    }

private:
    std::optional<Function> g_left_;
    std::optional<Function> g_right_;

    static void enforce_node(std::size_t node,
                             double value,
                             linalg::SparseMatrix& A,
                             linalg::Vector& b);
};



} // namespace boundary
} // namespace fem