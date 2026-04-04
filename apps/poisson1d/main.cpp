#include <fem/fem.hpp>
#include <cmath>

static constexpr double PI = 3.141592653589793238462643383279502884;

int main()
{
    const std::size_t n_elements = 80;
    const double x_min = 0.0;
    const double x_max = 1.0;
    const double h = (x_max - x_min) / n_elements;

    std::vector<double> nodes;
    nodes.reserve(n_elements + 1);
    for (std::size_t i = 0; i <= n_elements; ++i)
        nodes.push_back(x_min + i * h);

    auto mesh = fem::core::Mesh(nodes);

    fem::problems::Poisson1D problem_sin{
        [](double x){ return PI * PI * std::sin(PI * x); }
    };
    auto sin_bc = fem::boundary::DirichletBC(
        [](double) { return 0.0; },
        [](double) { return 0.0; }
    );

    fem::problems::Poisson1D problem_linear{
        [](double) { return 0.0; }
    };
    auto linear_bc = fem::boundary::DirichletBC(
        [](double) { return 0.0; },
        [](double) { return 1.0; }
    );

    auto problem = problem_linear;
    auto bc      = linear_bc;

    const auto fe   = fem::discretization::element::LagrangeP1_1D();
    const auto quad = fem::discretization::quadrature::GaussLegendre1D(2);

    const auto res = fem::Driver::solve(mesh, problem, fe, quad, {&bc});

    fem::io::write_results(mesh, res, problem, fe, quad, {&bc});

    return 0;
}