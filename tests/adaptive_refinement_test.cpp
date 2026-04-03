#include "tests/test_assert.hpp"
#include "tests/test_helpers.hpp"

#include "fem/solve/driver.hpp"
#include "fem/problems/problem.hpp"
#include "fem/boundary/dirichlet_bc.hpp"
#include "fem/discretization/element/lagrange_p1_1d.hpp"
#include "fem/discretization/quadrature/gauss_legendre_1d.hpp"
#include "fem/analysis/jump_indicator.hpp"
#include "fem/analysis/refine.hpp"

#include <cmath>
#include <cstddef>

namespace {

// u(x) = sin(pi*x), stark variierend in der Mitte
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
    fem::discretization::element::LagrangeP1_1D fe;
    fem::discretization::quadrature::GaussLegendre1D quad(3);

    fem::boundary::DirichletBC bc(
        [](double) { return 0.0; },
        [](double) { return 0.0; }
    );

    const auto u_exact = [](double x) { return std::sin(M_PI * x); };

    // Grobes Ausgangsmesh
    auto mesh = fem::tests::make_uniform_mesh(0.0, 1.0, 8);

    const auto u0 = fem::Driver::solve(mesh, problem, fe, quad, {&bc});
    const double e0 = fem::tests::l2_error(mesh, u0, fe, quad, u_exact);

    // Ein Verfeinerungsschritt
    const auto est = fem::analysis::jump_indicator(mesh, u0, problem);

    fem::test::require(est.indicators.size() == mesh.n_elements(),
        "indicators size matches n_elements");
    fem::test::require(est.global_error > 0.0,
        "global error > 0");

    // Alle Indikatoren nicht-negativ
    for (std::size_t e = 0; e < est.indicators.size(); ++e) {
        fem::test::require(est.indicators[e] >= 0.0,
            "indicator[" + std::to_string(e) + "] >= 0");
    }

    const auto mesh_refined = fem::analysis::refine(mesh, est, 0.5);

    // Verfeinertes Mesh hat mehr Elemente
    fem::test::require(mesh_refined.n_elements() > mesh.n_elements(),
        "refined mesh has more elements");

    // L2-Fehler nimmt ab nach Verfeinerung
    const auto u1 = fem::Driver::solve(mesh_refined, problem, fe, quad, {&bc});
    const double e1 = fem::tests::l2_error(mesh_refined, u1, fe, quad, u_exact);

    fem::test::require(e1 < e0,
        "L2 error decreases after refinement");

    // Mehrere Verfeinerungsschritte konvergieren
    auto mesh_curr = mesh_refined;
    double e_curr  = e1;

    for (std::size_t step = 0; step < 4; ++step) {
        const auto u_curr = fem::Driver::solve(mesh_curr, problem, fe, quad, {&bc});
        const auto est_curr = fem::analysis::jump_indicator(mesh_curr, u_curr, problem);
        const auto mesh_next = fem::analysis::refine(mesh_curr, est_curr, 0.5);
        const auto u_next = fem::Driver::solve(mesh_next, problem, fe, quad, {&bc});
        const double e_next = fem::tests::l2_error(mesh_next, u_next, fe, quad, u_exact);

        fem::test::require(e_next < e_curr,
            "error decreases at step " + std::to_string(step));

        mesh_curr = mesh_next;
        e_curr    = e_next;
    }

    return 0;
}