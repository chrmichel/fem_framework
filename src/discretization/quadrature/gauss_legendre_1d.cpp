#include "fem/discretization/quadrature/gauss_legendre_1d.hpp"

#include <stdexcept>
#include <sstream>

namespace fem::discretization::quadrature {

GaussLegendre1D::GaussLegendre1D(std::size_t n_points) {
    init(n_points);
}

double GaussLegendre1D::point(std::size_t q) const {
    if (q >= m_points.size()) {
        std::ostringstream oss;
        oss << "GaussLegendre1D::point: index q=" << q
            << " out of range (size=" << m_points.size() << ")";
        throw std::out_of_range(oss.str());
    }
    return m_points[q];
}

double GaussLegendre1D::weight(std::size_t q) const {
    if (q >= m_weights.size()) {
        std::ostringstream oss;
        oss << "GaussLegendre1D::weight: index q=" << q
            << " out of range (size=" << m_weights.size() << ")";
        throw std::out_of_range(oss.str());
    }
    return m_weights[q];
}

std::string GaussLegendre1D::name() const {
    std::ostringstream oss;
    oss << "gauss_legendre_1d(n=" << m_n_points << ")";
    return oss.str();
}

void GaussLegendre1D::init(std::size_t n_points) {
    m_n_points = n_points;
    m_points.clear();
    m_weights.clear();

    // Rules on [-1,1]. Values are standard Gauss-Legendre nodes/weights.
    // Phase-2 scope: support small n (1..5).
    switch (n_points) {
    case 1: {
        m_points  = { 0.0 };
        m_weights = { 2.0 };
        break;
    }
    case 2: {
        const double a = 0.5773502691896257645091487805019574556476; // 1/sqrt(3)
        m_points  = { -a, +a };
        m_weights = { 1.0, 1.0 };
        break;
    }
    case 3: {
        const double a = 0.7745966692414833770358530799564799221666;
        m_points  = { -a, 0.0, +a };
        m_weights = { 0.5555555555555555555555555555555555555556,
                      0.8888888888888888888888888888888888888889,
                      0.5555555555555555555555555555555555555556 };
        break;
    }
    case 4: {
        const double a = 0.3399810435848562648026657591032446872006;
        const double b = 0.8611363115940525752239464888928095050957;
        m_points  = { -b, -a, +a, +b };
        m_weights = { 0.3478548451374538573730639492219994072353,
                      0.6521451548625461426269360507780005927647,
                      0.6521451548625461426269360507780005927647,
                      0.3478548451374538573730639492219994072353 };
        break;
    }
    case 5: {
        const double a = 0.5384693101056830910363144207002088049673;
        const double b = 0.9061798459386639927976268782993929651257;
        m_points  = { -b, -a, 0.0, +a, +b };
        m_weights = { 0.2369268850561890875142640407199173626433,
                      0.4786286704993664680412915148356381929123,
                      0.5688888888888888888888888888888888888889,
                      0.4786286704993664680412915148356381929123,
                      0.2369268850561890875142640407199173626433 };
        break;
    }
    default: {
        std::ostringstream oss;
        oss << "GaussLegendre1D: unsupported n_points=" << n_points
            << " (supported: 1..5 in Phase 2)";
        throw std::invalid_argument(oss.str());
    }
    }

    if (m_points.size() != m_weights.size()) {
        throw std::runtime_error("GaussLegendre1D internal error: points/weights size mismatch");
    }
}

} // namespace fem::discretization::quadrature