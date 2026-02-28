#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <fem/fem.hpp>

static constexpr double PI = 3.141592653589793238462643383279502884;

static double u_exact(double x)
{
    return std::sin(PI * x);
}

static double f_rhs(double x)
{
    return PI * PI * std::sin(PI * x);
}

static double max_node_error(const fem::Mesh& mesh, const fem::linalg::Vector& u)
{
    double max_err = 0.0;
    for (std::size_t i = 0; i < mesh.n_nodes(); ++i) {
        const double x = mesh.node(i);
        max_err = std::max(max_err, std::abs(u(i) - u_exact(x)));
    }
    return max_err;
}

int main()
{
    fem::solve::Driver driver;

    std::vector<std::size_t> Ns = {10, 20, 40, 80};
    std::vector<double> errs;

    for (auto N : Ns) {
        fem::problems::Poisson1D problem(
            0.0, 1.0, N,
            [](double x){ return f_rhs(x); },  // rhs
            [](double){ return 1.0; },         // diffusion
            [](double){ return 0.0; },         // g_left
            [](double){ return 0.0; }          // g_right
        );

        auto u = driver.solve(problem);
        const auto& mesh = problem.mesh();

        const double e = max_node_error(mesh, u);
        errs.push_back(e);

        std::cout << "N=" << N << " max_node_error=" << e << "\n";
    }

    // Fehler muss sinken
    for (std::size_t k = 1; k < errs.size(); ++k) {
        assert(errs[k] < errs[k-1]);
    }

    // Grobe Konvergenzrate checken: e(h)/e(h/2) sollte ~4 sein (tolerant)
    // Wir verlangen mindestens Faktor 2.5 (robuster gegen Rand-/Norm-Effekte).
    for (std::size_t k = 1; k < errs.size(); ++k) {
        const double ratio = errs[k-1] / errs[k];
        std::cout << "ratio=" << ratio << "\n";
        assert(ratio > 2.5);
    }

    std::cout << "PoissonConvergenceTest passed.\n";
    return 0;
}