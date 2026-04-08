#include "tests/test_assert.hpp"
#include "tests/test_helpers.hpp"

#include "fem/core/mesh_utils.hpp"
#include "fem/core/mesh_impl.hpp"
#include "fem/assembly/assembler.hpp"
#include "fem/assembly/sparsity.hpp"
#include "fem/discretization/element/lagrange_p1_2d.hpp"
#include "fem/discretization/quadrature/triangle_quadrature.hpp"
#include "fem/problems/problem.hpp"
#include "fem/linalg/vector.hpp"

#include <cmath>

namespace {

class Poisson2D final : public fem::problems::Problem {
public:
    static constexpr double PI = 3.141592653589793238462643383279502884;

    double rhs(double x, double y) const override {
        return 2.0 * PI * PI * std::sin(PI * x) * std::sin(PI * y);
    }
    double exact(double x, double y) const override {
        return std::sin(PI * x) * std::sin(PI * y);
    }
    std::string name() const override { return "Poisson2D"; }
};

} // namespace

int main() {
    const auto mesh = fem::core::create_triangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 4, 4);

    Poisson2D problem;
    fem::discretization::element::LagrangeP1_2D fe;
    fem::discretization::quadrature::TriangleQuadrature quad(3);

    // Sparsity-Pattern aufbauen
    auto A = fem::assembly::build_sparse_matrix(mesh, fe);
    fem::linalg::Vector b(mesh.n_nodes());

    // Assembler aufrufen
    fem::assembly::Assembler<2> assembler(mesh, problem, fe, quad);
    assembler.assemble(A, b);

    // Matrix ist nicht leer
    fem::test::require(A.nnz() > 0, "A has non-zero entries");

    // Diagonale ist positiv (SPD vor BC)
    for (std::size_t i = 0; i < mesh.n_nodes(); ++i) {
        fem::test::require(A(i, i) > 0.0,
            "diagonal entry A(" + std::to_string(i) + "," +
            std::to_string(i) + ") > 0");
    }

    // RHS ist nicht überall null
    double b_norm = 0.0;
    for (std::size_t i = 0; i < mesh.n_nodes(); ++i)
        b_norm += b(i) * b(i);
    fem::test::require(b_norm > 0.0, "RHS is not zero");

    // Zeilensumme der Steifigkeitsmatrix ist nahe null
    // (Eigenschaft der FEM-Steifigkeitsmatrix für Laplace)
    for (std::size_t i = 0; i < mesh.n_nodes(); ++i) {
        double row_sum = 0.0;
        for (std::size_t k = A.row_ptr()[i]; k < A.row_ptr()[i+1]; ++k)
            row_sum += A.values()[k];
        fem::test::require_close(row_sum, 0.0, 1e-10,
            "row sum of stiffness matrix near zero for row " +
            std::to_string(i));
    }

    return 0;
}