#include "tests/test_assert.hpp"
#include "tests/test_helpers.hpp"

#include "fem/solve/driver.hpp"
#include "fem/io/output.hpp"
#include "fem/problems/problem.hpp"
#include "fem/boundary/dirichlet_bc.hpp"
#include "fem/discretization/element/lagrange_p1_1d.hpp"
#include "fem/discretization/quadrature/gauss_legendre_1d.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace fs = std::filesystem;

namespace {

class OutputTestProblem final : public fem::problems::Problem {
public:
    double rhs(double) const override {
        return 0.0;
    }

    double exact(double x) const override {
        return x;
    }

    std::string name() const override {
        return "OutputTestProblem";
    }

    fem::io::Meta meta() const override {
        return {
            {"case", std::string("linear_dirichlet")},
            {"alpha", 1.0},
            {"order", std::int64_t{1}},
            {"stable", true}
        };
    }
};

std::string read_file(const fs::path& path) {
    std::ifstream in(path);
    if (!in) {
        fem::test::fail("could not open file: " + path.string());
    }

    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

bool contains(const std::string& text, const std::string& needle) {
    return text.find(needle) != std::string::npos;
}

} // namespace

int main() {
    const auto mesh = fem::tests::make_uniform_mesh(0.0, 1.0, 10);

    OutputTestProblem problem;
    fem::discretization::element::LagrangeP1_1D fe;
    fem::discretization::quadrature::GaussLegendre1D quad(2);

    fem::boundary::DirichletBC<1> bc(
        [](auto) { return 0.0; },
        [](auto) { return 1.0; }
    );

    const auto u = fem::Driver<1>::solve(mesh, problem, fe, quad, {&bc});

    const fs::path base_dir = "test_results_output";
    if (fs::exists(base_dir)) {
        fs::remove_all(base_dir);
    }

    const fs::path run_dir =
        fem::io::write_results(mesh, u, problem, fe, quad, {&bc}, base_dir.string());

    fem::test::require(fs::exists(run_dir), "result directory exists");
    fem::test::require(fs::is_directory(run_dir), "result path is a directory");

    const fs::path solution_file = run_dir / "solution.csv";
    const fs::path meta_file     = run_dir / "meta.json";

    fem::test::require(fs::exists(solution_file), "solution.csv exists");
    fem::test::require(fs::exists(meta_file), "meta.json exists");

    const std::string solution_text = read_file(solution_file);
    const std::string meta_text     = read_file(meta_file);

    // CSV sanity checks
    fem::test::require(contains(solution_text, "node_id,x,u"), "solution.csv has header");
    fem::test::require(contains(solution_text, "0,0"), "solution.csv contains first node row");

    // JSON/meta sanity checks
    fem::test::require(contains(meta_text, "\"problem\""), "meta.json contains problem key");
    fem::test::require(contains(meta_text, "\"problem_meta\""), "meta.json contains problem_meta key");
    fem::test::require(contains(meta_text, "\"finite_element\""), "meta.json contains finite_element key");
    fem::test::require(contains(meta_text, "\"quadrature\""), "meta.json contains quadrature key");
    fem::test::require(contains(meta_text, "\"boundary_condition\""), "meta.json contains boundary_condition key");
    fem::test::require(contains(meta_text, "\"n_nodes\""), "meta.json contains n_nodes key");
    fem::test::require(contains(meta_text, "\"n_elements\""), "meta.json contains n_elements key");

    fem::test::require(contains(meta_text, "OutputTestProblem"), "meta.json contains problem name");
    fem::test::require(contains(meta_text, "lagrange_p1_1d"), "meta.json contains finite element name");
    fem::test::require(contains(meta_text, "gauss_legendre_1d"), "meta.json contains quadrature name");
    fem::test::require(contains(meta_text, "DirichletBC"), "meta.json contains BC name");

    // problem_meta contents
    fem::test::require(contains(meta_text, "\"case\""), "meta.json contains meta key 'case'");
    fem::test::require(contains(meta_text, "\"linear_dirichlet\""), "meta.json contains meta value 'linear_dirichlet'");
    fem::test::require(contains(meta_text, "\"alpha\""), "meta.json contains meta key 'alpha'");
    fem::test::require(contains(meta_text, "\"order\""), "meta.json contains meta key 'order'");
    fem::test::require(contains(meta_text, "\"stable\""), "meta.json contains meta key 'stable'");

    fs::remove_all(base_dir);

    return 0;
}