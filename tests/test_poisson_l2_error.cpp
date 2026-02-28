#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

#include <fem/fem.hpp>

static constexpr double PI = 3.141592653589793238462643383279502884;

// exakte Lösung
static double u_exact(double x)
{
    return std::sin(PI * x);
}

// RHS: -u'' = pi^2 sin(pi x)
static double f_rhs(double x)
{
    return PI * PI * std::sin(PI * x);
}

// P1-Shape-Funktionen auf Referenzelement [-1,1]
static double phi(int a, double xi)
{
    if (a == 0) return 0.5 * (1.0 - xi);
    if (a == 1) return 0.5 * (1.0 + xi);
    return 0.0;
}

// berechnet L2-Fehler mit 2-Punkt Gauss pro Element
static double l2_error(const fem::Mesh& mesh, const fem::linalg::Vector& u)
{
    // 2-Punkt Gauss auf [-1,1]
    const double xi[2] = { -0.5773502691896257, 0.5773502691896257 };
    const double w[2]  = {  1.0,                 1.0 };

    double integral = 0.0;

    for (std::size_t e = 0; e < mesh.n_elements(); ++e) {
        auto conn = mesh.element(e);
        const std::size_t i0 = conn[0];
        const std::size_t i1 = conn[1];

        const double x0 = mesh.node(i0);
        const double x1 = mesh.node(i1);
        const double J  = (x1 - x0) / 2.0;

        for (int q = 0; q < 2; ++q) {
            const double x = 0.5 * (x0 + x1) + J * xi[q];

            // u_h(x) = sum_a u_i_a * N_a(xi)
            const double uh =
                u(i0) * phi(0, xi[q]) +
                u(i1) * phi(1, xi[q]);

            const double diff = u_exact(x) - uh;
            integral += w[q] * J * diff * diff;
        }
    }

    return std::sqrt(integral);
}

int main()
{
    fem::solve::Driver driver;

    std::vector<std::size_t> Ns = {10, 20, 40, 80};
    std::vector<double> errs;
    errs.reserve(Ns.size());

    for (auto N : Ns) {
        fem::problems::Poisson1D problem(
            0.0, 1.0, N,
            [](double x) { return f_rhs(x); },  // rhs
            [](double) { return 1.0; },         // diffusion
            [](double) { return 0.0; },         // g_left
            [](double) { return 0.0; }          // g_right
        );

        auto u = driver.solve(problem);
        const auto& mesh = problem.mesh();

        const double e = l2_error(mesh, u);
        errs.push_back(e);

        std::cout << "N=" << N << " L2_error=" << e << "\n";
    }

    // Fehler muss sinken
    for (std::size_t k = 1; k < errs.size(); ++k) {
        assert(errs[k] < errs[k - 1]);
    }

    // Für P1 erwartet man im L2: O(h^2) => Ratio ~ 4 beim Halbieren von h
    // Tolerant prüfen: Ratio > 3 ist robust.
    for (std::size_t k = 1; k < errs.size(); ++k) {
        const double ratio = errs[k - 1] / errs[k];
        std::cout << "ratio=" << ratio << "\n";
        assert(ratio > 3.0);
    }

    std::cout << "PoissonL2ErrorTest passed.\n";
    return 0;
}