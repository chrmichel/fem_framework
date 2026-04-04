#include "fem/linalg/solver.hpp"
#include "fem/linalg/checks.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace fem::linalg {

// solver.cpp — neue Überladung für SparseMatrix
linalg::Vector solve(const linalg::SparseMatrix& A,
                     const linalg::Vector& b)
{
    const std::size_t n = b.size();
    linalg::Vector x(n);  // Startwert: 0

    // Matrix-Vektor-Produkt
    auto matvec = [&](const linalg::Vector& v) {
        linalg::Vector result(n);
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t k = A.row_ptr()[i]; k < A.row_ptr()[i+1]; ++k)
                sum += A.values()[k] * v(A.col_idx()[k]);
            result(i) = sum;
        }
        return result;
    };

    // Skalarprodukt
    auto dot = [&](const linalg::Vector& u, const linalg::Vector& v) {
        double sum = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            sum += u(i) * v(i);
        return sum;
    };

    linalg::Vector r(b);       // r = b - A*x = b (x=0)
    linalg::Vector p(r);       // p = r
    double rr = dot(r, r);

    for (std::size_t iter = 0; iter < 10 * n; ++iter) {
        if (rr < 1e-24) break;

        linalg::Vector Ap = matvec(p);
        const double pAp   = dot(p, Ap);
        const double alpha = rr / pAp;

        for (std::size_t i = 0; i < n; ++i) {
            x(i) += alpha * p(i);
            r(i) -= alpha * Ap(i);
        }

        const double rr_new = dot(r, r);
        const double beta   = rr_new / rr;
        rr = rr_new;

        for (std::size_t i = 0; i < n; ++i)
            p(i) = r(i) + beta * p(i);
    }

    return x;
}

} // namespace fem::linalg