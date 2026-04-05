#pragma once
#include "fem/discretization/element/finite_element.hpp"
#include <cstddef>
#include <string>

namespace fem::discretization::element {

class LagrangeP1_2D final : public FiniteElement {
public:
    std::size_t num_dofs() const override { return 3; }

    // 1D-Methoden werfen — dieses Element ist 2D
    double shape(std::size_t i, double xi) const override {
        throw std::logic_error("LagrangeP1_2D: use shape(i, xi, eta)");
    }
    double dshape_dxi(std::size_t i, double xi) const override {
        throw std::logic_error("LagrangeP1_2D: use dshape_dxi(i, xi, eta)");
    }

    double shape(std::size_t i, double xi, double eta) const override;
    double dshape_dxi(std::size_t i, double xi, double eta) const override;
    double dshape_deta(std::size_t i, double xi, double eta) const override;

    std::string name() const override { return "lagrange_p1_2d"; }
};

} // namespace fem::discretization::element