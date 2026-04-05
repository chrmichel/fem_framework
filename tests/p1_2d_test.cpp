#include "tests/test_assert.hpp"
#include "fem/discretization/element/lagrange_p1_2d.hpp"

#include <cmath>

int main() {
    fem::discretization::element::LagrangeP1_2D fe;

    fem::test::require(fe.num_dofs() == 3, "num_dofs == 3");

    // Interpolationseigenschaft: phi_i(xi_j, eta_j) == delta_ij
    const double xi[3]  = {0.0, 1.0, 0.0};
    const double eta[3] = {0.0, 0.0, 1.0};

    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            const double expected = (i == j) ? 1.0 : 0.0;
            fem::test::require_close(
                fe.shape(i, xi[j], eta[j]), expected, 1e-14,
                "shape(" + std::to_string(i) + ", xi[" + std::to_string(j) + "]) == delta");
        }
    }

    // Partition of Unity: sum phi_i == 1
    for (auto [x, e] : std::initializer_list<std::pair<double,double>>{
            {0.1, 0.2}, {0.3, 0.3}, {0.5, 0.1}, {0.0, 0.0}, {0.25, 0.25}}) {
        double sum = 0.0;
        for (std::size_t i = 0; i < 3; ++i)
            sum += fe.shape(i, x, e);
        fem::test::require_close(sum, 1.0, 1e-14,
            "partition of unity at (" + std::to_string(x) + "," + std::to_string(e) + ")");
    }

    // Ableitungen: konstant, analytisch prüfen
    for (auto [x, e] : std::initializer_list<std::pair<double,double>>{
            {0.1, 0.2}, {0.3, 0.1}}) {

        // dshape_dxi
        fem::test::require_close(fe.dshape_dxi(0, x, e), -1.0, 1e-14, "dphi0/dxi == -1");
        fem::test::require_close(fe.dshape_dxi(1, x, e),  1.0, 1e-14, "dphi1/dxi ==  1");
        fem::test::require_close(fe.dshape_dxi(2, x, e),  0.0, 1e-14, "dphi2/dxi ==  0");

        // dshape_deta
        fem::test::require_close(fe.dshape_deta(0, x, e), -1.0, 1e-14, "dphi0/deta == -1");
        fem::test::require_close(fe.dshape_deta(1, x, e),  0.0, 1e-14, "dphi1/deta ==  0");
        fem::test::require_close(fe.dshape_deta(2, x, e),  1.0, 1e-14, "dphi2/deta ==  1");
    }

    // 1D-Methoden werfen
    bool threw = false;
    try { fe.shape(0, 0.5); } catch (const std::logic_error&) { threw = true; }
    fem::test::require(threw, "shape(i, xi) throws for 2D element");

    threw = false;
    try { fe.dshape_dxi(0, 0.5); } catch (const std::logic_error&) { threw = true; }
    fem::test::require(threw, "dshape_dxi(i, xi) throws for 2D element");

    return 0;
}