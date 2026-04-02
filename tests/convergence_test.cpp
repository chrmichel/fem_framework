#include "tests/test_assert.hpp"
#include "tests/test_helpers.hpp"

#include "fem/solve/driver.hpp"
#include "fem/problems/problem.hpp"
#include "fem/boundary/dirichlet_bc.hpp"
#include "fem/discretization/element/lagrange_p1_1d.hpp"
#include "fem/discretization/quadrature/gauss_legendre_1d.hpp"

#include "fem/analysis/l2_error.hpp"

#include <cmath>
#include <cstddef>

namespace {

class SinProblem final : public fem::problems::Problem {
public:
    double rhs(double x) const override {
        return M_PI * M_PI * std::sin(M_PI * x);
    }

    double exact(double x) const override {
        return std::sin(M_PI * x);
    }

    std::string name() const override {
        return "SinProblem";
    }
};

} // namespace

int main() {
    SinProblem problem;
    fem::discretization::element::LagrangeP1_1D fe;
    fem::discretization::quadrature::GaussLegendre1D quad(3);

    fem::boundary::DirichletBC bc(
        [](double) { return 0.0; },
        [](double) { return 0.0; }
    );

    const auto u_exact = [](double x) { return std::sin(M_PI * x); };

    const std::size_t n1 = 20;
    const std::size_t n2 = 40;
    const std::size_t n3 = 80;

    const auto mesh1 = fem::tests::make_uniform_mesh(0.0, 1.0, n1);
    const auto mesh2 = fem::tests::make_uniform_mesh(0.0, 1.0, n2);
    const auto mesh3 = fem::tests::make_uniform_mesh(0.0, 1.0, n3);

    const auto u1 = fem::Driver::solve(mesh1, problem, fe, quad, {&bc});
    const auto u2 = fem::Driver::solve(mesh2, problem, fe, quad, {&bc});
    const auto u3 = fem::Driver::solve(mesh3, problem, fe, quad, {&bc});

    const double e1 = fem::analysis::l2_error(mesh1, u1, fe, quad, u_exact);
    const double e2 = fem::analysis::l2_error(mesh2, u2, fe, quad, u_exact);
    const double e3 = fem::analysis::l2_error(mesh3, u3, fe, quad, u_exact);

    const double h1 = 1.0 / static_cast<double>(n1);
    const double h2 = 1.0 / static_cast<double>(n2);
    const double h3 = 1.0 / static_cast<double>(n3);

    const double p12 = std::log(e1 / e2) / std::log(h1 / h2);
    const double p23 = std::log(e2 / e3) / std::log(h2 / h3);

    fem::test::require(p12 > 1.5, "first observed convergence rate > 1.5");
    fem::test::require(p23 > 1.5, "second observed convergence rate > 1.5");

    return 0;
}