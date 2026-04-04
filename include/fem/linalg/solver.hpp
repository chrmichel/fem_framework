#pragma once

#include "vector.hpp"
#include "sparse_matrix.hpp"

namespace fem::linalg {

// Phase 3: Sparse solver (Conjugate Gradient)
Vector solve(const SparseMatrix& A, const Vector& b);

} // namespace fem::linalg