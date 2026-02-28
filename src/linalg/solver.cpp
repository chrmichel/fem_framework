#include "fem/linalg/solver.hpp"
#include "fem/linalg/checks.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace fem::linalg {

Vector solve(const Matrix& A_in, const Vector& b_in)
{
    require_square(A_in);
    require_same_size(A_in, b_in);

    const std::size_t n = A_in.rows();

    // Make local copies (solver modifies them)
    Matrix A = A_in;
    Vector b = b_in;

    // Gaussian elimination with partial pivoting
    for (std::size_t k = 0; k < n; ++k) {
        // Pivot search
        std::size_t piv = k;
        double max_abs = std::abs(A(k, k));
        for (std::size_t i = k + 1; i < n; ++i) {
            const double v = std::abs(A(i, k));
            if (v > max_abs) {
                max_abs = v;
                piv = i;
            }
        }

        if (max_abs == 0.0) {
            throw std::runtime_error("solve: singular matrix (zero pivot).");
        }

        // Swap rows k <-> piv in A and b
        if (piv != k) {
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(A(k, j), A(piv, j));
            }
            std::swap(b(k), b(piv));
        }

        // Elimination below pivot
        const double Akk = A(k, k);
        for (std::size_t i = k + 1; i < n; ++i) {
            const double factor = A(i, k) / Akk;

            // Set explicitly to 0 to reduce numeric noise
            A(i, k) = 0.0;

            for (std::size_t j = k + 1; j < n; ++j) {
                A(i, j) -= factor * A(k, j);
            }
            b(i) -= factor * b(k);
        }
    }

    // Back substitution
    Vector x(n, 0.0);
    for (std::size_t ii = 0; ii < n; ++ii) {
        const std::size_t i = n - 1 - ii;

        double sum = b(i);
        for (std::size_t j = i + 1; j < n; ++j) {
            sum -= A(i, j) * x(j);
        }

        const double Aii = A(i, i);
        if (Aii == 0.0) {
            throw std::runtime_error("solve: singular matrix (zero diagonal in back-substitution).");
        }

        x(i) = sum / Aii;
    }

    return x;
}

} // namespace fem::linalg