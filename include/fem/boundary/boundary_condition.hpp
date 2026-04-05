#pragma once
#include <string>

namespace fem {
namespace linalg { class Vector; class SparseMatrix; }
namespace core { template<int Dim> class Mesh; }

namespace boundary {

template<int Dim>
class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;

    virtual void apply(linalg::SparseMatrix& A,
                       linalg::Vector& b,
                       const core::Mesh<Dim>& mesh) const = 0;

    virtual std::string name() const { return "BoundaryCondition"; }
};

} // namespace boundary
} // namespace fem