#pragma once
#include <string>

namespace fem {
namespace linalg { class Vector; class SparseMatrix; }
namespace core { class Mesh; }

namespace boundary {

class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;

    virtual void apply(linalg::SparseMatrix& A,
                       linalg::Vector& b,
                       const core::Mesh& mesh) const = 0;

    virtual std::string name() const { return "BoundaryCondition"; }
};

} // namespace boundary
} // namespace fem