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
    for (std::size_t e = 0; e < m_nodes.size() - 1; ++e) {
        m_connectivity.push_back({e, e + 1});
    }
}

Mesh::Mesh(std::vector<double> nodes,
           std::vector<std::vector<std::size_t>> connectivity)
    : m_nodes(std::move(nodes)),
      m_connectivity(std::move(connectivity))
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
    // validate connectivity
    for (std::size_t e = 0; e < m_connectivity.size(); ++e) {
        if (m_connectivity[e].size() < 2) {
            std::ostringstream oss;
            oss << "Mesh: element " << e << " has fewer than 2 nodes.";
            throw std::invalid_argument(oss.str());
        }
        for (std::size_t local = 0; local < m_connectivity[e].size(); ++local) {
            if (m_connectivity[e][local] >= m_nodes.size()) {
                std::ostringstream oss;
                oss << "Mesh: element " << e << ", local node " << local
                    << " has index " << m_connectivity[e][local]
                    << " out of range (n_nodes=" << m_nodes.size() << ").";
                throw std::out_of_range(oss.str());
            }
        }
    }
}

std::size_t Mesh::n_nodes() const {
    return m_nodes.size();
}

std::size_t Mesh::n_elements() const {
    return m_connectivity.size();
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

    if (local_index >= m_connectivity[e].size()) {
        std::ostringstream oss;
        oss << "Mesh::element_node: local_index "
            << local_index << " invalid for element " << e << ", (ndofs="
            << m_connectivity[e].size() << ")";
        throw std::out_of_range(oss.str());
    }

    return m_connectivity[e][local_index];
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