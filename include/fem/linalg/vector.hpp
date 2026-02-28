#pragma once

#include <cstddef>
#include <vector>
#include <initializer_list>

namespace fem::linalg {

class Vector {
public:
    using Index = std::size_t;

    Vector() = default;
    explicit Vector(Index n, double value = 0.0);
    Vector(std::initializer_list<double> init);

    Index size() const noexcept;

    // Elementzugriff (mit Bounds-Check)
    double& operator()(Index i);
    double  operator()(Index i) const;

    // raw access (für Solver/IO, optional)
    double* data() noexcept;
    const double* data() const noexcept;

    void set_zero() noexcept;
    void fill(double value) noexcept;

private:
    std::vector<double> data_;

    void check_index(Index i) const;
};

} // namespace fem::linalg