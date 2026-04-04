#pragma once
#include <string>

namespace fem {
namespace linalg { class Vector; class SparseMatrix; }
namespace core { template<int Dim> class Mesh; }

namespace boundary {

class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;

    virtual void apply(linalg::SparseMatrix& A,
                       linalg::Vector& b,
                       const core::Mesh<1>& mesh) const = 0;

    virtual std::string name() const { return "BoundaryCondition"; }
};

} // namespace boundary
} // namespace fem