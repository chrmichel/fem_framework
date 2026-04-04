#include "tests/test_assert.hpp"
#include "tests/test_helpers.hpp"

#include "fem/assembly/sparsity.hpp"
#include "fem/discretization/element/lagrange_p1_1d.hpp"
#include "fem/discretization/element/lagrange_p2_1d.hpp"
#include "fem/core/mesh_utils.hpp"

int main() {
    fem::discretization::element::LagrangeP1_1D fe_p1;
    fem::discretization::element::LagrangeP2_1D fe_p2;

    // P1: 3 Knoten, 2 Elemente -> Tridiagonal -> 7 Einträge
    const auto mesh_p1 = fem::tests::make_uniform_mesh(0.0, 1.0, 2);
    const auto A_p1 = fem::assembly::build_sparse_matrix(mesh_p1, fe_p1);

    fem::test::require(A_p1.rows() == 3, "P1: rows == 3");
    fem::test::require(A_p1.nnz()  == 7, "P1: nnz == 7");

    // P1: Diagonale und Nebendiagonalen existieren
    // (0,0), (0,1), (1,0), (1,1), (1,2), (2,1), (2,2)
    bool threw = false;
    try { (void)A_p1(0, 2); } catch (const std::out_of_range&) { threw = true; }
    fem::test::require(threw, "P1: (0,2) not in pattern");

    threw = false;
    try { (void)A_p1(0, 0); } catch (...) { threw = true; }
    fem::test::require(!threw, "P1: (0,0) in pattern");

    // P2: 5 Knoten, 2 Elemente -> 17 Einträge
    const auto mesh_p2 = fem::core::create_p2_mesh_1d(0.0, 1.0, 2);
    const auto A_p2 = fem::assembly::build_sparse_matrix(mesh_p2, fe_p2);

    fem::test::require(A_p2.rows() == 5, "P2: rows == 5");
    fem::test::require(A_p2.nnz()  == 17, "P2: nnz == 17");

    return 0;
}