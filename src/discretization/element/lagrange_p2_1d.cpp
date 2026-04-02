#include "fem/discretization/element/lagrange_p2_1d.hpp"

#include <stdexcept>
#include <sstream>

namespace fem::discretization::element {

double LagrangeP2_1D::shape(std::size_t i, double xi) const {
    switch (i) {
    case 0:
        return 0.5 * xi * (xi - 1.0);
    case 1:
        return (1.0 + xi) * (1.0 - xi);
    case 2:
        return 0.5 * xi * (1.0 + xi);
    default: {
        std::ostringstream oss;
        oss << "LagrangeP2_1D::shape: local index i=" << i
            << " out of range (num_dofs=3)";
        throw std::out_of_range(oss.str());
    }
    }
}

double LagrangeP2_1D::dshape_dxi(std::size_t i, double xi) const {
    switch (i) {
    case 0:
        return xi - 0.5;
    case 1:
        return -2.0 * xi;
    case 2:
        return xi + 0.5;
    default: {
        std::ostringstream oss;
        oss << "LagrangeP2_1D::dshape_dxi: local index i=" << i
            << " out of range (num_dofs=3)";
        throw std::out_of_range(oss.str());
    }
    }
}

} // namespace fem::discretization::element