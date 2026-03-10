#include "fem/discretization/element/lagrange_p1_1d.hpp"

#include <stdexcept>
#include <sstream>

namespace fem::discretization::element {

double LagrangeP1_1D::shape(std::size_t i, double xi) const {
    switch (i) {
    case 0:
        return 0.5 * (1.0 - xi);
    case 1:
        return 0.5 * (1.0 + xi);
    default: {
        std::ostringstream oss;
        oss << "LagrangeP1_1D::shape: local index i=" << i
            << " out of range (num_dofs=2)";
        throw std::out_of_range(oss.str());
    }
    }
}

double LagrangeP1_1D::dshape_dxi(std::size_t i, double /*xi*/) const {
    switch (i) {
    case 0:
        return -0.5;
    case 1:
        return +0.5;
    default: {
        std::ostringstream oss;
        oss << "LagrangeP1_1D::dshape_dxi: local index i=" << i
            << " out of range (num_dofs=2)";
        throw std::out_of_range(oss.str());
    }
    }
}

} // namespace fem::discretization::element