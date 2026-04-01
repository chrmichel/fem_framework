#include <fem/fem.hpp>
#include <cmath>
#include <iostream>

static constexpr double PI = 3.141592653589793238462643383279502884;

int main()
{
    auto nodes = std::vector<double>();
    const std::size_t n_elements = 80;
    const double x_min = 0.0;
    const double x_max = 1.0;
    const double h = (x_max - x_min) / n_elements;
    for (std::size_t i = 0; i <= n_elements; ++i) {
        nodes.push_back(x_min + i * h);
    }
    auto mesh = fem::core::Mesh(nodes);
    // Problem definieren + Labels (damit meta.json sinnvoll ist)
    fem::problems::Poisson1D problem_sin(
        [](double x){ return PI * PI * std::sin(PI * x); } // rhs
        //"pi^2*sin(pi*x)", "1", "0", "0"
    );
    // auto u_exact_sin = [](double x){ return std::sin(PI * x); };  // optional exact
    auto sin_bc = fem::boundary::DirichletBC(
        [](double x){ return 0.0; }, // g_left
        [](double x){ return 0.0; }  // g_right
    );

    fem::problems::Poisson1D problem_linear(
        [](double x){ return 0.0; }                          // rhs
    );
    // auto u_exact_linear = [](double x){ return x; };  
    auto linear_bc = fem::boundary::DirichletBC(
        [](double x){ return 0.0; }, // g_left
        [](double x){ return 1.0; }  // g_right
    );

    auto problem = problem_linear; // switch problem here
    // auto u_exact = u_exact_linear; // optional exact
    auto bc = linear_bc; // switch BC here

    auto res = fem::Driver::solve(
        mesh, problem,
        fem::discretization::element::LagrangeP1_1D(),
        fem::discretization::quadrature::GaussLegendre1D(2),
        bc
    );

    // fully automatic meta.json from problem + mesh + framework defaults
    fem::io::write_results(
        mesh, res, problem, 
        fem::discretization::element::LagrangeP1_1D(),
        fem::discretization::quadrature::GaussLegendre1D(2),
        bc
    );

    //std::cout << "Wrote results to: " << out_dir << "\n";
    return 0;
}