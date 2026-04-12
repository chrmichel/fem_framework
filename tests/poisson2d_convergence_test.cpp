#include "tests/test_assert.hpp"
#include "tests/test_helpers.hpp"

#include "fem/core/mesh_utils.hpp"
#include "fem/core/mesh_impl.hpp"
#include "fem/solve/driver.hpp"
#include "fem/problems/problem.hpp"
#include "fem/boundary/dirichlet_bc.hpp"
#include "fem/discretization/element/lagrange_p1_2d.hpp"
#include "fem/discretization/quadrature/triangle_quadrature.hpp"

#include <cmath>
#include <cstddef>

namespace {

class SinSinProblem final : public fem::problems::Problem {
public:
    static constexpr double PI = 3.141592653589793238462643383279502884;

    double rhs(double x, double y) const override {
        return 2.0 * PI * PI * std::sin(PI * x) * std::sin(PI * y);
    }
    double exact(double x, double y) const override {
        return std::sin(PI * x) * std::sin(PI * y);
    }
    std::string name() const override { return "SinSinProblem"; }
};

} // namespace

int main() {
    SinSinProblem problem;
    fem::discretization::element::LagrangeP1_2D fe;
    fem::discretization::quadrature::TriangleQuadrature quad(3);

    // Dirichlet u=0 auf allen Rändern
    fem::boundary::DirichletBC<2> bc(
        [](auto p) { return 0.0; },
        std::nullopt
    );

    const auto u_exact = [](auto p) {
        constexpr double PI = 3.141592653589793238462643383279502884;
        return std::sin(PI * p[0]) * std::sin(PI * p[1]);
    };

    // Drei Gitter: 4x4, 8x8, 16x16
    const auto mesh1 = fem::core::create_triangle_mesh_2d(0.0,1.0, 0.0,1.0, 4,  4);
    const auto mesh2 = fem::core::create_triangle_mesh_2d(0.0,1.0, 0.0,1.0, 8,  8);
    const auto mesh3 = fem::core::create_triangle_mesh_2d(0.0,1.0, 0.0,1.0, 16, 16);

    const auto u1 = fem::Driver<2>::solve(mesh1, problem, fe, quad, {&bc});
    const auto u2 = fem::Driver<2>::solve(mesh2, problem, fe, quad, {&bc});
    const auto u3 = fem::Driver<2>::solve(mesh3, problem, fe, quad, {&bc});

    const double e1 = fem::tests::l2_error<2>(mesh1, u1, fe, quad, u_exact);
    const double e2 = fem::tests::l2_error<2>(mesh2, u2, fe, quad, u_exact);
    const double e3 = fem::tests::l2_error<2>(mesh3, u3, fe, quad, u_exact);

    const double h1 = 1.0 / 4.0;
    const double h2 = 1.0 / 8.0;
    const double h3 = 1.0 / 16.0;

    const double p12 = std::log(e1 / e2) / std::log(h1 / h2);
    const double p23 = std::log(e2 / e3) / std::log(h2 / h3);

    fem::test::require(e2 < e1, "error decreases n1->n2");
    fem::test::require(e3 < e2, "error decreases n2->n3");
    fem::test::require(p12 > 1.5, "convergence rate p12 > 1.5");
    fem::test::require(p23 > 1.5, "convergence rate p23 > 1.5");

    return 0;
}