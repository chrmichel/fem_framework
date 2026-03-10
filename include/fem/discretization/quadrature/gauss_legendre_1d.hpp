#pragma once

#include "fem/discretization/quadrature/quadrature_rule.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace fem::discretization::quadrature {

/**
 * Gauss-Legendre quadrature on [-1,1].
 *
 * Phase 2 scope:
 *   Provide n-point rules for small n (e.g. 1..5).
 */
class GaussLegendre1D final : public QuadratureRule {
public:
    explicit GaussLegendre1D(std::size_t n_points);

    std::size_t size() const override { return m_points.size(); }
    double point(std::size_t q) const override;
    double weight(std::size_t q) const override;

    std::string name() const override;

private:
    std::vector<double> m_points;
    std::vector<double> m_weights;
    std::size_t m_n_points{0};

    void init(std::size_t n_points);
};

} // namespace fem::discretization::quadrature