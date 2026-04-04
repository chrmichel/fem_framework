#include "tests/test_assert.hpp"
#include "fem/linalg/sparse_matrix.hpp"

#include <stdexcept>

int main() {
    // Einfache 3x3 Tridiagonalmatrix
    // [ 2 -1  0 ]
    // [-1  2 -1 ]
    // [ 0 -1  2 ]
    std::vector<std::size_t> row_ptr = {0, 2, 5, 7};
    std::vector<std::size_t> col_idx = {0, 1, 0, 1, 2, 1, 2};

    fem::linalg::SparseMatrix A(3, row_ptr, col_idx);

    fem::test::require(A.rows() == 3, "rows == 3");
    fem::test::require(A.cols() == 3, "cols == 3");
    fem::test::require(A.nnz()  == 7, "nnz == 7");

    // Einträge setzen
    A(0, 0) =  2.0; A(0, 1) = -1.0;
    A(1, 0) = -1.0; A(1, 1) =  2.0; A(1, 2) = -1.0;
    A(2, 1) = -1.0; A(2, 2) =  2.0;

    fem::test::require_close(A(0, 0),  2.0, 1e-14, "A(0,0) == 2");
    fem::test::require_close(A(0, 1), -1.0, 1e-14, "A(0,1) == -1");
    fem::test::require_close(A(1, 1),  2.0, 1e-14, "A(1,1) == 2");
    fem::test::require_close(A(2, 2),  2.0, 1e-14, "A(2,2) == 2");

    // set_zero
    A.set_zero();
    fem::test::require_close(A(0, 0), 0.0, 1e-14, "after set_zero: A(0,0) == 0");
    fem::test::require_close(A(1, 2), 0.0, 1e-14, "after set_zero: A(1,2) == 0");

    // Zugriff auf nicht-existenten Eintrag wirft
    bool threw = false;
    try {
        double v = A(0, 2);
        (void)v;
    } catch (const std::out_of_range&) {
        threw = true;
    }
    fem::test::require(threw, "access to non-pattern entry throws");

    return 0;
}