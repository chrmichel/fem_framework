#include "fem/core/mesh.hpp"

#include <stdexcept>
#include <sstream>

namespace fem::core {

Mesh::Mesh(std::vector<double> nodes)
    : m_nodes(std::move(nodes))
{
    if (m_nodes.size() < 2) {
        throw std::invalid_argument("Mesh: need at least 2 nodes.");
    }
    for (std::size_t i = 1; i < m_nodes.size(); ++i) {
        if (m_nodes[i] <= m_nodes[i-1]) {
            std::ostringstream oss;
            oss << "Mesh: nodes must be strictly increasing, but node "
                << i << " has coordinate " << m_nodes[i]
                << " <= previous node " << m_nodes[i-1];
            throw std::invalid_argument(oss.str());
        }
    }
}

std::size_t Mesh::n_nodes() const {
    return m_nodes.size();
}

std::size_t Mesh::n_elements() const {
    return m_nodes.size() - 1;
}

double Mesh::node(std::size_t i) const {
    if (i >= m_nodes.size()) {
        std::ostringstream oss;
        oss << "Mesh::node: index " << i
            << " out of range (num_nodes=" << m_nodes.size() << ")";
        throw std::out_of_range(oss.str());
    }
    return m_nodes[i];
}

std::size_t Mesh::element_node(std::size_t e,
                               std::size_t local_index) const
{
    if (e >= n_elements()) {
        std::ostringstream oss;
        oss << "Mesh::element_node: element index " << e
            << " out of range (n_elements=" << n_elements() << ")";
        throw std::out_of_range(oss.str());
    }

    if (local_index > 1) {
        std::ostringstream oss;
        oss << "Mesh::element_node: local_index "
            << local_index << " invalid for P1 element (must be 0 or 1)";
        throw std::out_of_range(oss.str());
    }

    return e + local_index;
}
// --- NEW ---

std::size_t Mesh::left_boundary_node() const {
    return 0;
}

std::size_t Mesh::right_boundary_node() const {
    return n_nodes() - 1;
}

std::array<std::size_t, 2> Mesh::boundary_nodes() const {
    return { left_boundary_node(), right_boundary_node() };
}

double Mesh::x_left() const {
    return node(left_boundary_node());
}

double Mesh::x_right() const {
    return node(right_boundary_node());
}
} // namespace fem::core