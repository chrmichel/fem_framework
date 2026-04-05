#include "fem/discretization/quadrature/triangle_quadrature.hpp"
#include <stdexcept>
#include <sstream>

namespace fem::discretization::quadrature {

TriangleQuadrature::TriangleQuadrature(std::size_t n_points) {
    init(n_points);
}

std::array<double, 2> TriangleQuadrature::point2d(std::size_t q) const {
    if (q >= m_points.size()) {
        std::ostringstream oss;
        oss << "TriangleQuadrature::point2d: index q=" << q
            << " out of range (size=" << m_points.size() << ")";
        throw std::out_of_range(oss.str());
    }
    return m_points[q];
}

double TriangleQuadrature::weight(std::size_t q) const {
    if (q >= m_weights.size()) {
        std::ostringstream oss;
        oss << "TriangleQuadrature::weight: index q=" << q
            << " out of range (size=" << m_weights.size() << ")";
        throw std::out_of_range(oss.str());
    }
    return m_weights[q];
}

std::string TriangleQuadrature::name() const {
    std::ostringstream oss;
    oss << "triangle_quadrature(n=" << m_n_points << ")";
    return oss.str();
}

void TriangleQuadrature::init(std::size_t n_points) {
    m_n_points = n_points;
    m_points.clear();
    m_weights.clear();

    switch (n_points) {
    case 1: {
        // Schwerpunktregel — exakt für lineare Polynome
        m_points  = { {1.0/3.0, 1.0/3.0} };
        m_weights = { 0.5 };
        break;
    }
    case 3: {
        // 3-Punkt-Regel — exakt für quadratische Polynome
        m_points  = { {1.0/6.0, 1.0/6.0},
                      {2.0/3.0, 1.0/6.0},
                      {1.0/6.0, 2.0/3.0} };
        m_weights = { 1.0/6.0, 1.0/6.0, 1.0/6.0 };
        break;
    }
    default: {
        std::ostringstream oss;
        oss << "TriangleQuadrature: unsupported n_points=" << n_points
            << " (supported: 1, 3)";
        throw std::invalid_argument(oss.str());
    }
    }
}

} // namespace fem::discretization::quadrature