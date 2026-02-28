#include "fem/linalg/matrix.hpp"
#include <stdexcept>

namespace fem::linalg {

Matrix::Matrix(Index rows, Index cols, double value)
    : rows_(rows), cols_(cols), data_(rows * cols, value)
{
}

Matrix::Index Matrix::rows() const noexcept { return rows_; }
Matrix::Index Matrix::cols() const noexcept { return cols_; }

double& Matrix::operator()(Index i, Index j)
{
    check_index(i, j);
    return data_[i * cols_ + j];
}

double Matrix::operator()(Index i, Index j) const
{
    check_index(i, j);
    return data_[i * cols_ + j];
}

double* Matrix::data() noexcept { return data_.data(); }
const double* Matrix::data() const noexcept { return data_.data(); }

void Matrix::set_zero() noexcept
{
    fill(0.0);
}

void Matrix::fill(double value) noexcept
{
    for (auto& v : data_) v = value;
}

void Matrix::set_row_zero(Index i)
{
    check_row(i);
    for (Index j = 0; j < cols_; ++j) {
        data_[i * cols_ + j] = 0.0;
    }
}

void Matrix::set_col_zero(Index j)
{
    check_col(j);
    for (Index i = 0; i < rows_; ++i) {
        data_[i * cols_ + j] = 0.0;
    }
}

void Matrix::check_index(Index i, Index j) const
{
    if (i >= rows_ || j >= cols_) {
        throw std::out_of_range("Matrix index out of range.");
    }
}

void Matrix::check_row(Index i) const
{
    if (i >= rows_) {
        throw std::out_of_range("Matrix row index out of range.");
    }
}

void Matrix::check_col(Index j) const
{
    if (j >= cols_) {
        throw std::out_of_range("Matrix col index out of range.");
    }
}

} // namespace fem::linalg