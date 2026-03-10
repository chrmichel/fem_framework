#pragma once

#include "boundary_condition.hpp"
#include <functional>
#include <cstddef>

namespace fem {
class Mesh;

namespace boundary {

class DirichletBC final : public BoundaryCondition {
public:
    using Function = std::function<double(double)>;

    DirichletBC(Function g_left, Function g_right);

    void apply(linalg::Matrix& A,
               linalg::Vector& b,
               const core::Mesh& mesh) const override;

private:
    Function g_left_;
    Function g_right_;

    static void enforce_node(std::size_t node,
                             double value,
                             linalg::Matrix& A,
                             linalg::Vector& b);
    
    virtual std::string name() const override {
        return "DirichletBC";
    }
};



} // namespace boundary
} // namespace fem