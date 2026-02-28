#pragma once

#include <cstddef>
#include <vector>

namespace fem::linalg {

class Matrix {
public:
    using Index = std::size_t;

    Matrix() = default;
    Matrix(Index rows, Index cols, double value = 0.0);

    Index rows() const noexcept;
    Index cols() const noexcept;

    // Elementzugriff (mit Bounds-Check)
    double& operator()(Index i, Index j);
    double  operator()(Index i, Index j) const;

    // raw access (row-major)
    double* data() noexcept;
    const double* data() const noexcept;

    void set_zero() noexcept;
    void fill(double value) noexcept;

    // Hilfen für BC/Elimination
    void set_row_zero(Index i);
    void set_col_zero(Index j);

private:
    Index rows_{0};
    Index cols_{0};
    std::vector<double> data_; // row-major: data_[i*cols_ + j]

    void check_index(Index i, Index j) const;
    void check_row(Index i) const;
    void check_col(Index j) const;
};

} // namespace fem::linalg