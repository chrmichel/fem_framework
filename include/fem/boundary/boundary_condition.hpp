#pragma once

namespace fem {
namespace linalg { class Matrix; class Vector; }
class Mesh;

namespace boundary {

class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;

    virtual void apply(linalg::Matrix& A,
                       linalg::Vector& b,
                       const Mesh& mesh) const = 0;
};

} // namespace boundary
} // namespace fem