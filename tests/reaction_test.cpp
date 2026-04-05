#include "tests/test_assert.hpp"
#include "tests/test_helpers.hpp"

#include "fem/solve/driver.hpp"
#include "fem/problems/problem.hpp"
#include "fem/boundary/dirichlet_bc.hpp"
#include "fem/discretization/element/lagrange_p1_1d.hpp"
#include "fem/discretization/quadrature/gauss_legendre_1d.hpp"

#include <cmath>
#include <cstddef>

namespace {

class NoReactionProblem final : public fem::problems::Problem {
public:
    double rhs(double) const override {
        return 1.0;
    }

    double reaction(double) const override {
        return 0.0;
    }

    std::string name() const override {
        return "NoReactionProblem";
    }
};

class WithReactionProblem final : public fem::problems::Problem {
public:
    double rhs(double) const override {
        return 1.0;
    }

    double reaction(double) const override {
        return 10.0;
    }

    std::string name() const override {
        return "WithReactionProblem";
    }
};

} // namespace

int main() {
    const auto mesh = fem::tests::make_uniform_mesh(0.0, 1.0, 40);

    fem::discretization::element::LagrangeP1_1D fe;
    fem::discretization::quadrature::GaussLegendre1D quad(3);

    fem::boundary::DirichletBC<1> bc(
        [](double) { return 0.0; },
        [](double) { return 0.0; }
    );

    NoReactionProblem p0;
    WithReactionProblem p1;

    const auto u0 = fem::Driver<1>::solve(mesh, p0, fe, quad, {&bc});
    const auto u1 = fem::Driver<1>::solve(mesh, p1, fe, quad, {&bc});

    double diff_norm = 0.0;
    for (std::size_t i = 0; i < mesh.n_nodes(); ++i) {
        const double d = u0(i) - u1(i);
        diff_norm += d * d;
    }

    diff_norm = std::sqrt(diff_norm);

    fem::test::require(diff_norm > 1e-8, "reaction term changes the solution");

    return 0;
}