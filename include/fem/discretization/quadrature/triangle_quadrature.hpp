#pragma once
#include "fem/discretization/quadrature/quadrature_rule.hpp"
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace fem::discretization::quadrature {

/**
 * Quadratur auf dem Referenzdreieck {(0,0),(1,0),(0,1)}.
 *
 * Exakt für Polynome bis Grad n_points-1.
 * Phase 4 scope: n_points = 1, 3.
 *
 * Integral: integral_{Referenzdreieck} f dA
 *         = sum_q f(xi_q, eta_q) * w_q
 *
 * Fläche des Referenzdreiecks = 0.5,
 * Gewichte summieren sich zu 0.5.
 */
class TriangleQuadrature final : public QuadratureRule {
public:
    explicit TriangleQuadrature(std::size_t n_points);

    std::size_t size() const override { return m_points.size(); }

    // 1D-Methode wirft
    double point(std::size_t q) const override {
        throw std::logic_error("TriangleQuadrature: use point2d(q)");
    }

    std::array<double, 2> point2d(std::size_t q) const override;
    double weight(std::size_t q) const override;
    std::string name() const override;

private:
    std::vector<std::array<double, 2>> m_points;
    std::vector<double>                m_weights;
    std::size_t                        m_n_points{0};

    void init(std::size_t n_points);
};

} // namespace fem::discretization::quadrature