#pragma once

#include <cstddef>
#include <string>

namespace fem::discretization::quadrature {

/**
 * Quadrature rule on a 1D reference segment.
 *
 * Convention for Phase 2:
 *   Reference segment: xi in [-1, 1]
 *
 * Provides quadrature points xi_q and weights w_q such that
 *   ∫_{-1}^{1} g(xi) dxi ≈ Σ_q g(xi_q) * w_q
 */
class QuadratureRule {
public:
    virtual ~QuadratureRule() = default;

    /// Number of quadrature points
    virtual std::size_t size() const = 0;

    /// q-th quadrature point in reference coordinate xi (xi in [-1,1])
    virtual double point(std::size_t q) const = 0;

    /// q-th quadrature weight
    virtual double weight(std::size_t q) const = 0;

    /// Human-readable identifier (useful for logging/tests)
    virtual std::string name() const = 0;
};

} // namespace fem::discretization::quadrature