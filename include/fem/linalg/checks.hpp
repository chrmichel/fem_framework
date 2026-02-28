#pragma once
#include "matrix.hpp"
#include "vector.hpp"
#include <stdexcept>

namespace fem::linalg {

inline void require_square(const Matrix& A)
{
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Matrix must be square.");
    }
}

inline void require_same_size(const Matrix& A, const Vector& b)
{
    if (A.rows() != b.size()) {
        throw std::invalid_argument("Size mismatch: A.rows() != b.size().");
    }
}

} // namespace fem::linalg