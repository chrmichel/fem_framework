#include "tests/test_assert.hpp"
#include "fem/discretization/quadrature/triangle_quadrature.hpp"

#include <cmath>

int main() {
    // --- 1-Punkt-Regel ---
    fem::discretization::quadrature::TriangleQuadrature quad1(1);

    fem::test::require(quad1.size() == 1, "1-point: size == 1");

    // Gewichte summieren sich zu 0.5 (Fläche des Referenzdreiecks)
    double sum_w1 = 0.0;
    for (std::size_t q = 0; q < quad1.size(); ++q)
        sum_w1 += quad1.weight(q);
    fem::test::require_close(sum_w1, 0.5, 1e-14, "1-point: weights sum to 0.5");

    // Punkte liegen im Referenzdreieck (xi >= 0, eta >= 0, xi+eta <= 1)
    for (std::size_t q = 0; q < quad1.size(); ++q) {
        const auto p = quad1.point2d(q);
        fem::test::require(p[0] >= 0.0, "1-point: xi >= 0");
        fem::test::require(p[1] >= 0.0, "1-point: eta >= 0");
        fem::test::require(p[0] + p[1] <= 1.0 + 1e-14, "1-point: xi+eta <= 1");
    }

    // Exakt für Konstante: integral_{T} 1 dA = 0.5
    {
        double result = 0.0;
        for (std::size_t q = 0; q < quad1.size(); ++q)
            result += 1.0 * quad1.weight(q);
        fem::test::require_close(result, 0.5, 1e-14, "1-point: integrates 1 exactly");
    }

    // Exakt für lineare Funktion: integral_{T} xi dA = 1/6
    {
        double result = 0.0;
        for (std::size_t q = 0; q < quad1.size(); ++q)
            result += quad1.point2d(q)[0] * quad1.weight(q);
        fem::test::require_close(result, 1.0/6.0, 1e-14,
            "1-point: integrates xi exactly");
    }

    // --- 3-Punkt-Regel ---
    fem::discretization::quadrature::TriangleQuadrature quad3(3);

    fem::test::require(quad3.size() == 3, "3-point: size == 3");

    // Gewichte summieren sich zu 0.5
    double sum_w3 = 0.0;
    for (std::size_t q = 0; q < quad3.size(); ++q)
        sum_w3 += quad3.weight(q);
    fem::test::require_close(sum_w3, 0.5, 1e-14, "3-point: weights sum to 0.5");

    // Punkte liegen im Referenzdreieck
    for (std::size_t q = 0; q < quad3.size(); ++q) {
        const auto p = quad3.point2d(q);
        fem::test::require(p[0] >= 0.0, "3-point: xi >= 0");
        fem::test::require(p[1] >= 0.0, "3-point: eta >= 0");
        fem::test::require(p[0] + p[1] <= 1.0 + 1e-14, "3-point: xi+eta <= 1");
    }

    // Exakt für Konstante
    {
        double result = 0.0;
        for (std::size_t q = 0; q < quad3.size(); ++q)
            result += 1.0 * quad3.weight(q);
        fem::test::require_close(result, 0.5, 1e-14, "3-point: integrates 1 exactly");
    }

    // Exakt für xi
    {
        double result = 0.0;
        for (std::size_t q = 0; q < quad3.size(); ++q)
            result += quad3.point2d(q)[0] * quad3.weight(q);
        fem::test::require_close(result, 1.0/6.0, 1e-14,
            "3-point: integrates xi exactly");
    }

    // Exakt für xi*eta: integral_{T} xi*eta dA = 1/24
    {
        double result = 0.0;
        for (std::size_t q = 0; q < quad3.size(); ++q) {
            const auto p = quad3.point2d(q);
            result += p[0] * p[1] * quad3.weight(q);
        }
        fem::test::require_close(result, 1.0/24.0, 1e-14,
            "3-point: integrates xi*eta exactly");
    }

    // 1D-Methode wirft
    bool threw = false;
    try { quad3.point(0); } catch (const std::logic_error&) { threw = true; }
    fem::test::require(threw, "point(q) throws for triangle quadrature");

    // Ungültige Punktanzahl wirft
    threw = false;
    try {
        fem::discretization::quadrature::TriangleQuadrature bad(2);
    } catch (const std::invalid_argument&) { threw = true; }
    fem::test::require(threw, "unsupported n_points throws");

    return 0;
}