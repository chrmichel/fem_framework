// include/fem/core/mesh_impl.hpp
#pragma once
#include "mesh.hpp"
#include <stdexcept>
#include <sstream>

namespace fem::core {

template<int Dim>
Mesh<Dim>::Mesh(std::vector<Point> nodes)
    : m_nodes(std::move(nodes))
{
    validate_nodes();
    for (std::size_t e = 0; e < m_nodes.size() - 1; ++e)
        m_connectivity.push_back({e, e + 1});
    m_left_boundary_node  = 0;
    m_right_boundary_node = m_nodes.size() - 1;
}

template<int Dim>
Mesh<Dim>::Mesh(std::vector<Point> nodes,
                std::vector<std::vector<std::size_t>> connectivity)
    : m_nodes(std::move(nodes)),
      m_connectivity(std::move(connectivity))
{
    validate_nodes();
    validate_connectivity();
    m_left_boundary_node  = 0;
    m_right_boundary_node = m_nodes.size() - 1;
}

template<int Dim>
std::size_t Mesh<Dim>::n_nodes() const {
    return m_nodes.size();
}

template<int Dim>
std::size_t Mesh<Dim>::n_elements() const {
    return m_connectivity.size();
}

template<int Dim>
std::size_t Mesh<Dim>::nodes_per_element(std::size_t e) const {
    return m_connectivity[e].size();
}

template<int Dim>
typename Mesh<Dim>::Point Mesh<Dim>::node(std::size_t i) const {
    if (i >= m_nodes.size()) {
        std::ostringstream oss;
        oss << "Mesh::node: index " << i
            << " out of range (n_nodes=" << m_nodes.size() << ")";
        throw std::out_of_range(oss.str());
    }
    return m_nodes[i];
}

template<int Dim>
std::size_t Mesh<Dim>::element_node(std::size_t e,
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
        oss << "Mesh::element_node: local_index " << local_index
            << " invalid for element " << e
            << " (ndofs=" << m_connectivity[e].size() << ")";
        throw std::out_of_range(oss.str());
    }
    return m_connectivity[e][local_index];
}

template<int Dim>
std::size_t Mesh<Dim>::left_boundary_node() const {
    return m_left_boundary_node;
}

template<int Dim>
std::size_t Mesh<Dim>::right_boundary_node() const {
    return m_right_boundary_node;
}

template<int Dim>
std::array<std::size_t, 2> Mesh<Dim>::boundary_nodes() const {
    return { m_left_boundary_node, m_right_boundary_node };
}

template<int Dim>
void Mesh<Dim>::set_boundary_nodes(std::size_t left, std::size_t right) {
    if (left >= m_nodes.size())
        throw std::out_of_range("Mesh::set_boundary_nodes: left index out of range");
    if (right >= m_nodes.size())
        throw std::out_of_range("Mesh::set_boundary_nodes: right index out of range");
    m_left_boundary_node  = left;
    m_right_boundary_node = right;
}

template<int Dim>
void Mesh<Dim>::validate_nodes() const {
    if (m_nodes.size() < 2)
        throw std::invalid_argument("Mesh: need at least 2 nodes.");
}

template<int Dim>
void Mesh<Dim>::validate_connectivity() const {
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

} // namespace fem::core