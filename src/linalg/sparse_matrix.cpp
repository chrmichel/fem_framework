#include "fem/linalg/sparse_matrix.hpp"

#include <algorithm>
#include <stdexcept>
#include <sstream>

namespace fem::linalg {

SparseMatrix::SparseMatrix(std::size_t n_rows,
                           std::vector<std::size_t> row_ptr,
                           std::vector<std::size_t> col_idx)
    : m_rows(n_rows),
      m_row_ptr(std::move(row_ptr)),
      m_col_idx(std::move(col_idx)),
      m_values(m_col_idx.size(), 0.0)
{
    if (m_row_ptr.size() != m_rows + 1) {
        throw std::invalid_argument(
            "SparseMatrix: row_ptr must have n_rows+1 entries");
    }
}

std::size_t SparseMatrix::rows() const { return m_rows; }
std::size_t SparseMatrix::cols() const { return m_rows; }
std::size_t SparseMatrix::nnz()  const { return m_values.size(); }

void SparseMatrix::set_zero() {
    std::fill(m_values.begin(), m_values.end(), 0.0);
}

const std::vector<std::size_t>& SparseMatrix::row_ptr() const { return m_row_ptr; }
const std::vector<std::size_t>& SparseMatrix::col_idx() const { return m_col_idx; }
const std::vector<double>&      SparseMatrix::values()  const { return m_values; }

std::size_t SparseMatrix::find(std::size_t i, std::size_t j) const {
    const std::size_t start = m_row_ptr[i];
    const std::size_t end   = m_row_ptr[i + 1];

    for (std::size_t k = start; k < end; ++k) {
        if (m_col_idx[k] == j) return k;
    }

    std::ostringstream oss;
    oss << "SparseMatrix: entry (" << i << "," << j << ") not in sparsity pattern";
    throw std::out_of_range(oss.str());
}

double& SparseMatrix::operator()(std::size_t i, std::size_t j) {
    return m_values[find(i, j)];
}

double SparseMatrix::operator()(std::size_t i, std::size_t j) const {
    return m_values[find(i, j)];
}

} // namespace fem::linalg