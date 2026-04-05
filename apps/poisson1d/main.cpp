#include <fem/fem.hpp>
#include <cmath>

static constexpr double PI = 3.141592653589793238462643383279502884;

int main()
{
    const std::size_t n_elements = 80;
    const double x_min = 0.0;
    const double x_max = 1.0;

    auto mesh = fem::core::create_p1_mesh_1d(x_min, x_max, n_elements);

    fem::problems::Poisson1D problem_sin{
        [](double x){ return PI * PI * std::sin(PI * x); }
    };
    auto sin_bc = fem::boundary::DirichletBC<1>(
        [](double) { return 0.0; },
        [](double) { return 0.0; }
    );

    fem::problems::Poisson1D problem_linear{
        [](double) { return 0.0; }
    };
    auto linear_bc = fem::boundary::DirichletBC<1>(
        [](double) { return 0.0; },
        [](double) { return 1.0; }
    );

    auto problem = problem_linear;
    auto bc      = linear_bc;

    const auto fe   = fem::discretization::element::LagrangeP1_1D();
    const auto quad = fem::discretization::quadrature::GaussLegendre1D(2);

    const auto res = fem::Driver<1>::solve(mesh, problem, fe, quad, {&bc});

    fem::io::write_results(mesh, res, problem, fe, quad, {&bc});

    return 0;
}