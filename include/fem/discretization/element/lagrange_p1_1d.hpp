#pragma once

#include "fem/discretization/element/finite_element.hpp"

#include <cstddef>
#include <string>

namespace fem::discretization::element {

/**
 * 1D Lagrange P1 element on reference segment [-1,1].
 *
 * Local DoFs:
 *   i = 0  → xi = -1
 *   i = 1  → xi = +1
 *
 * Shape functions:
 *   phi0(xi) = (1 - xi)/2
 *   phi1(xi) = (1 + xi)/2
 *
 * Derivatives:
 *   dphi0/dxi = -1/2
 *   dphi1/dxi = +1/2
 */
class LagrangeP1_1D final : public FiniteElement {
public:
    std::size_t num_dofs() const override { return 2; }

    double shape(std::size_t i, double xi) const override;
    double dshape_dxi(std::size_t i, double xi) const override;

    std::string name() const override { return "lagrange_p1_1d"; }
};

} // namespace fem::discretization::element