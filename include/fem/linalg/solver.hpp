#pragma once

#include "matrix.hpp"
#include "vector.hpp"

namespace fem::linalg {

// Phase 1: Dense solver (Gaussian elimination with partial pivoting)
Vector solve(const Matrix& A, const Vector& b);

} // namespace fem::linalg