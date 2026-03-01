#include <fem/fem.hpp>
#include <cmath>
#include <iostream>

static constexpr double PI = 3.141592653589793238462643383279502884;

int main()
{
    // Problem definieren + Labels (damit meta.json sinnvoll ist)
    fem::problems::Poisson1D problem_sin(
        0.0, 1.0, 80,
        [](double x){ return PI * PI * std::sin(PI * x); }, // rhs
        [](double){ return 1.0; },                          // diffusion
        [](double){ return 0.0; },                          // g_left
        [](double){ return 0.0; },                          // g_right
        "pi^2*sin(pi*x)", "1", "0", "0"
    );
    auto u_exact_sin = [](double x){ return std::sin(PI * x); };  // optional exact
    fem::problems::Poisson1D problem_linear(
        0.0, 1.0, 80,
        [](double x){ return 0.0; },                          // rhs
        [](double){ return 1.0; },                          // diffusion
        [](double){ return 0.0; },                          // g_left
        [](double){ return 1.0; },                          // g_right
        "0", "1", "1", "0"
    );
    auto u_exact_linear = [](double x){ return x; };  

    auto problem = problem_linear; // switch problem here
    auto u_exact = u_exact_linear; // optional exact

    fem::solve::Driver driver;
    auto res = driver.solve_with_system(problem);

    const auto out_dir = fem::io::make_results_dir("results");

    fem::io::write_solution_csv(
        out_dir, problem.mesh(), res.u,
        u_exact
    );

    // fully automatic meta.json from problem + mesh + framework defaults
    fem::io::write_meta_json(out_dir, problem);

    std::cout << "Wrote results to: " << out_dir << "\n";
    return 0;
}