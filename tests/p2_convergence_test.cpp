#include "tests/test_assert.hpp"
#include "tests/test_helpers.hpp"

#include "fem/core/mesh_utils.hpp"
#include "fem/solve/driver.hpp"
#include "fem/problems/problem.hpp"
#include "fem/boundary/dirichlet_bc.hpp"
#include "fem/discretization/element/lagrange_p2_1d.hpp"
#include "fem/discretization/quadrature/gauss_legendre_1d.hpp"

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
    std::string name() const override { return "SinProblem"; }
};

} // namespace

int main() {
    SinProblem problem;
    fem::discretization::element::LagrangeP2_1D fe;
    fem::discretization::quadrature::GaussLegendre1D quad(3);

    fem::boundary::DirichletBC<1> bc(
        [](double) { return 0.0; },
        [](double) { return 0.0; }
    );

    const auto u_exact = [](double x) { return std::sin(M_PI * x); };

    const std::size_t n1 = 10;
    const std::size_t n2 = 20;
    const std::size_t n3 = 40;

    const auto mesh1 = fem::core::create_p2_mesh_1d(0.0, 1.0, n1);
    const auto mesh2 = fem::core::create_p2_mesh_1d(0.0, 1.0, n2);
    const auto mesh3 = fem::core::create_p2_mesh_1d(0.0, 1.0, n3);

    const auto u1 = fem::Driver<1>::solve(mesh1, problem, fe, quad, {&bc});
    const auto u2 = fem::Driver<1>::solve(mesh2, problem, fe, quad, {&bc});
    const auto u3 = fem::Driver<1>::solve(mesh3, problem, fe, quad, {&bc});

    const double e1 = fem::tests::l2_error(mesh1, u1, fe, quad, u_exact);
    const double e2 = fem::tests::l2_error(mesh2, u2, fe, quad, u_exact);
    const double e3 = fem::tests::l2_error(mesh3, u3, fe, quad, u_exact);

    const double h1 = 1.0 / static_cast<double>(n1);
    const double h2 = 1.0 / static_cast<double>(n2);
    const double h3 = 1.0 / static_cast<double>(n3);

    const double p12 = std::log(e1 / e2) / std::log(h1 / h2);
    const double p23 = std::log(e2 / e3) / std::log(h2 / h3);

    //debug print
    std::cerr << "e1=" << e1 << " e2=" << e2 << " e3=" << e3 << "\n";
    std::cerr << "p12=" << p12 << " p23=" << p23 << "\n";

    fem::test::require(e2 < e1, "error decreases with refinement (n1->n2)");
    fem::test::require(e3 < e2, "error decreases with refinement (n2->n3)");
    fem::test::require(p12 > 1.8, "convergence rate p12 > 1.8");
    fem::test::require(p23 > 1.8, "convergence rate p23 > 1.8");

    return 0;
}