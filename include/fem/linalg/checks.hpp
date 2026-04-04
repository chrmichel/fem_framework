#pragma once
#include "vector.hpp"
#include "sparse_matrix.hpp"
#include <stdexcept>

namespace fem::linalg {


inline void require_square(const linalg::SparseMatrix& A) {
    if (A.rows() != A.cols())
        throw std::invalid_argument("Matrix must be square");
}

inline void require_same_size(const linalg::SparseMatrix& A,
                              const linalg::Vector& b) {
    if (A.rows() != b.size())
        throw std::invalid_argument("Matrix and vector sizes do not match");
}

} // namespace fem::linalg