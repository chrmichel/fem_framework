#pragma once

#include <cstddef>
#include <string>

namespace fem::discretization::element {

/**
 * Finite element on a 1D reference segment.
 *
 * Reference segment convention (Phase 2):
 *   xi ∈ [-1, 1]
 *
 * IMPORTANT:
 *   This interface lives purely in reference space.
 *   No Mesh, no geometry mapping.
 */
class FiniteElement {
public:
    virtual ~FiniteElement() = default;

    /// Number of local degrees of freedom
    virtual std::size_t num_dofs() const = 0;

    /// Shape function phi_i(xi)
    virtual double shape(std::size_t i, double xi) const = 0;

    /// Reference derivative d/dxi phi_i(xi)
    virtual double dshape_dxi(std::size_t i, double xi) const = 0;

    /// Human-readable identifier
    virtual std::string name() const = 0;
};

} // namespace fem::discretization::element