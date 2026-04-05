#include "tests/test_assert.hpp"
#include "tests/test_helpers.hpp"

#include "fem/solve/driver.hpp"
#include "fem/discretization/element/lagrange_p1_1d.hpp"
#include "fem/discretization/quadrature/gauss_legendre_1d.hpp"

#include "fem/problems/problem.hpp"
#include "fem/boundary/dirichlet_bc.hpp"

#include "fem/analysis/l2_error.hpp"

#include <cmath>

namespace {

class SinPoisson1D final : public fem::problems::Problem {
public:
    double rhs(double x) const override { return M_PI * M_PI * std::sin(M_PI * x); }
    double diffusion(double) const override { return 1.0; }
};

} // namespace

int main() {
    SinPoisson1D problem;
    fem::discretization::element::LagrangeP1_1D fe;
    fem::discretization::quadrature::GaussLegendre1D quad(3);

    fem::boundary::DirichletBC<1> bc(
        [](double){ return 0.0; },
        [](double){ return 0.0; }
    );

    auto u_exact = [](double x) { return std::sin(M_PI * x); };

    const auto mesh1 = fem::tests::make_uniform_mesh(0.0, 1.0, 20);
    const auto mesh2 = fem::tests::make_uniform_mesh(0.0, 1.0, 40);

    const auto u1 = fem::Driver<1>::solve(mesh1, problem, fe, quad, {&bc});
    const auto u2 = fem::Driver<1>::solve(mesh2, problem, fe, quad, {&bc});

    const double e1 = fem::analysis::l2_error(mesh1, u1, fe, quad, u_exact);
    const double e2 = fem::analysis::l2_error(mesh2, u2, fe, quad, u_exact);

    fem::test::require(e1 > 0.0, "e1 > 0");
    fem::test::require(e2 > 0.0, "e2 > 0");
    fem::test::require(e2 < e1, "error decreases with refinement");

    return 0;
}