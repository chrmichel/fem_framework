#include "fem/linalg/vector.hpp"
#include <stdexcept>

namespace fem::linalg {

Vector::Vector(Index n, double value) : data_(n, value) {}

Vector::Vector(std::initializer_list<double> init) : data_(init) {}

Vector::Index Vector::size() const noexcept { return data_.size(); }

double& Vector::operator()(Index i)
{
    check_index(i);
    return data_[i];
}

double Vector::operator()(Index i) const
{
    check_index(i);
    return data_[i];
}

double* Vector::data() noexcept { return data_.data(); }
const double* Vector::data() const noexcept { return data_.data(); }

void Vector::set_zero() noexcept
{
    fill(0.0);
}

void Vector::fill(double value) noexcept
{
    for (auto& v : data_) v = value;
}

void Vector::check_index(Index i) const
{
    if (i >= data_.size()) {
        throw std::out_of_range("Vector index out of range.");
    }
}

} // namespace fem::linalg