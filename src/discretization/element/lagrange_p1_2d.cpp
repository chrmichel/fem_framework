#include "fem/discretization/element/lagrange_p1_2d.hpp"
#include <stdexcept>
#include <sstream>

namespace fem::discretization::element {

double LagrangeP1_2D::shape(std::size_t i, double xi, double eta) const {
    switch (i) {
    case 0: return 1.0 - xi - eta;
    case 1: return xi;
    case 2: return eta;
    default: {
        std::ostringstream oss;
        oss << "LagrangeP1_2D::shape: index " << i << " out of range";
        throw std::out_of_range(oss.str());
    }
    }
}

double LagrangeP1_2D::dshape_dxi(std::size_t i, double xi, double eta) const {
    switch (i) {
    case 0: return -1.0;
    case 1: return  1.0;
    case 2: return  0.0;
    default: {
        std::ostringstream oss;
        oss << "LagrangeP1_2D::dshape_dxi: index " << i << " out of range";
        throw std::out_of_range(oss.str());
    }
    }
}

double LagrangeP1_2D::dshape_deta(std::size_t i, double xi, double eta) const {
    switch (i) {
    case 0: return -1.0;
    case 1: return  0.0;
    case 2: return  1.0;
    default: {
        std::ostringstream oss;
        oss << "LagrangeP1_2D::dshape_deta: index " << i << " out of range";
        throw std::out_of_range(oss.str());
    }
    }
}

} // namespace fem::discretization::element