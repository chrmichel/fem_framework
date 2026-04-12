#include "tests/test_assert.hpp"
#include "tests/test_helpers.hpp"

#include "fem/solve/driver.hpp"
#include "fem/discretization/element/lagrange_p1_1d.hpp"
#include "fem/discretization/quadrature/gauss_legendre_1d.hpp"

#include "fem/problems/problem.hpp"
#include "fem/boundary/dirichlet_bc.hpp"

#include <cstddef>

namespace {

class ZeroRhsPoisson1D final : public fem::problems::Problem {
public:
    double rhs(double) const override { return 0.0; }
    double diffusion(double) const override { return 1.0; }
};

} // namespace

int main() {
    auto mesh = fem::tests::make_uniform_mesh(0.0, 1.0, 20);

    ZeroRhsPoisson1D problem;
    fem::discretization::element::LagrangeP1_1D fe;
    fem::discretization::quadrature::GaussLegendre1D quad(2);

    // u(0)=0, u(1)=1 => exact u(x)=x
    fem::boundary::DirichletBC<1> bc(
        [](auto){ return 0.0; },
        [](auto){ return 1.0; }
    );

    const auto u = fem::Driver<1>::solve(mesh, problem, fe, quad, {&bc});

    for (std::size_t i = 0; i < mesh.n_nodes(); ++i) {
        const double x = mesh.node(i)[0];
        fem::test::require_close(u(i), x, 1e-12, "nodal linear exactness u(i)=x");
    }

    return 0;
}