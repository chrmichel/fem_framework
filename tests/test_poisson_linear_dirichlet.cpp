#include <cassert>
#include <cmath>
#include <iostream>
#include <fem/fem.hpp>

static double linear_exact(double x, double gL, double gR)
{
    return gL + (gR - gL) * x;
}

int main()
{
    const double gL = 2.0;
    const double gR = -1.0;

    fem::problems::Poisson1D problem(
        0.0, 1.0, 25,
        [](double){ return 0.0; },  // rhs f(x)=0
        [](double){ return 1.0; },  // diffusion a(x)=1
        [=](double){ return gL; },  // left Dirichlet
        [=](double){ return gR; }   // right Dirichlet
    );

    fem::solve::Driver driver;
    auto u = driver.solve(problem);

    const auto& mesh = problem.mesh();

    // max-norm error
    double max_err = 0.0;
    for (std::size_t i = 0; i < mesh.n_nodes(); ++i) {
        const double x = mesh.node(i);
        const double ue = linear_exact(x, gL, gR);
        max_err = std::max(max_err, std::abs(u(i) - ue));
    }

    std::cout << "max_err = " << max_err << "\n";

    // sollte sehr klein sein (Pivoting + double => typ. ~1e-14..1e-12)
    assert(max_err < 1e-10);

    std::cout << "PoissonLinearDirichletTest passed.\n";
    return 0;
}