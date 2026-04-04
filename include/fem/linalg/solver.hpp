#pragma once

#include "matrix.hpp"
#include "vector.hpp"
#include "sparse_matrix.hpp"

namespace fem::linalg {

// Phase 1: Dense solver (Gaussian elimination with partial pivoting)
Vector solve(const Matrix& A, const Vector& b);
// Phase 3: Sparse solver (Conjugate Gradient)
Vector solve(const SparseMatrix& A, const Vector& b);

} // namespace fem::linalg