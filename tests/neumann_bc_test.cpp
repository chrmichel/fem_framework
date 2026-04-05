#include "tests/test_assert.hpp"
#include "tests/test_helpers.hpp"

#include "fem/solve/driver.hpp"
#include "fem/problems/problem.hpp"
#include "fem/boundary/dirichlet_bc.hpp"
#include "fem/boundary/neumann_bc.hpp"
#include "fem/discretization/element/lagrange_p1_1d.hpp"
#include "fem/discretization/quadrature/gauss_legendre_1d.hpp"

#include <cmath>
#include <cstddef>

namespace {

// u(x) = x^2
// -(u'')' = -2  =>  f(x) = -2
// u(0) = 0  (Dirichlet links)
// u'(1) = 2  (Neumann rechts)
class QuadraticProblem final : public fem::problems::Problem {
public:
    double rhs(double) const override { return -2.0; }
    double diffusion(double) const override { return 1.0; }
    std::string name() const override { return "QuadraticProblem"; }
};

} // namespace

int main() {
    const auto mesh = fem::tests::make_uniform_mesh(0.0, 1.0, 40);

    QuadraticProblem problem;
    fem::discretization::element::LagrangeP1_1D fe;
    fem::discretization::quadrature::GaussLegendre1D quad(2);

    fem::boundary::DirichletBC<1> dirichlet(
        [](double) { return 0.0; },  // u(0) = 0
        std::nullopt                  // rechter Rand: kein Dirichlet
    );
    fem::boundary::NeumannBC<1> neumann(
        std::nullopt,                 // linker Rand: kein Neumann
        [](double) { return 2.0; }   // u'(1) = 2
    );

    const auto u = fem::Driver<1>::solve(mesh, problem, fe, quad, {&dirichlet, &neumann});

    const auto u_exact = [](double x) { return x * x; };

    // L2-Konvergenz
    fem::discretization::quadrature::GaussLegendre1D quad_err(3);
    const double err = fem::tests::l2_error(mesh, u, fe, quad_err, u_exact);
    fem::test::require(err < 1e-3, "L2 error < 1e-3 for 40 elements");

    // Fehler nimmt mit Verfeinerung ab
    const auto mesh_fine = fem::tests::make_uniform_mesh(0.0, 1.0, 80);
    const auto u_fine = fem::Driver<1>::solve(mesh_fine, problem, fe, quad, {&neumann, &dirichlet});
    const double err_fine = fem::tests::l2_error(mesh_fine, u_fine, fe, quad_err, u_exact);
    fem::test::require(err_fine < err, "error decreases with refinement");

    return 0;
}