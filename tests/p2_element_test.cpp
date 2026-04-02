#include "tests/test_assert.hpp"
#include "fem/discretization/element/lagrange_p2_1d.hpp"

#include <cmath>

int main() {
    fem::discretization::element::LagrangeP2_1D fe;

    fem::test::require(fe.num_dofs() == 3, "num_dofs == 3");

    // Interpolationseigenschaft: phi_i(xi_j) == delta_ij
    const double xi[3] = {-1.0, 0.0, 1.0};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            const double expected = (i == j) ? 1.0 : 0.0;
            fem::test::require_close(fe.shape(i, xi[j]), expected, 1e-14,
                "shape(" + std::to_string(i) + ", xi[" + std::to_string(j) + "]) == delta");
        }
    }

    // Partition of Unity: sum phi_i(xi) == 1 für beliebige xi
    for (double x : {-0.7, -0.3, 0.0, 0.4, 0.9}) {
        double sum = 0.0;
        for (std::size_t i = 0; i < 3; ++i)
            sum += fe.shape(i, x);
        fem::test::require_close(sum, 1.0, 1e-14, "partition of unity at xi=" + std::to_string(x));
    }

    // Ableitungen numerisch verifizieren
    const double h = 1e-7;
    for (std::size_t i = 0; i < 3; ++i) {
        for (double x : {-0.7, -0.1, 0.0, 0.5, 0.8}) {
            const double fd = (fe.shape(i, x + h) - fe.shape(i, x - h)) / (2.0 * h);
            fem::test::require_close(fe.dshape_dxi(i, x), fd, 1e-6,
                "dshape_dxi(" + std::to_string(i) + ") numerisch bei xi=" + std::to_string(x));
        }
    }

    return 0;
}