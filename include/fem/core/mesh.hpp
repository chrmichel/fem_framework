#pragma once

#include <array>
#include <vector>
#include <cstddef>

namespace fem::core {

/**
 * 1D mesh for Phase 1/2.
 *
 * Nodes are stored in ascending order.
 * Connectivity (P1):
 *   element e connects nodes (e, e+1)
 */
class Mesh {
public:
    explicit Mesh(std::vector<double> nodes);
    Mesh(std::vector<double> nodes,
         std::vector<std::vector<std::size_t>> connectivity);

    std::size_t n_nodes() const;
    std::size_t n_elements() const;

    double node(std::size_t i) const;

    /// Returns global node index of local node in element e
    /// local_index must be 0 or 1 (P1 only)
    std::size_t element_node(std::size_t e,
                             std::size_t local_index) const;
 // --- NEW: boundary node ids (Phase 2: 1D) ---

    /// Global node id at left boundary (x = x_min)
    std::size_t left_boundary_node() const;

    /// Global node id at right boundary (x = x_max)
    std::size_t right_boundary_node() const;

    /// Convenience: both boundary node ids {left, right}
    std::array<std::size_t, 2> boundary_nodes() const;

    // Optional convenience: boundary coordinates
    double x_left() const;
    double x_right() const;
    std::size_t nodes_per_element(std::size_t e) const { return m_connectivity[e].size(); } 
    void set_boundary_nodes(std::size_t left_id, std::size_t right_id);
private:
    std::vector<double> m_nodes;
    std::vector<std::vector<std::size_t>> m_connectivity;
    std::size_t m_left_boundary_node{0}; // default to first node
    std::size_t m_right_boundary_node{0};
};

} // namespace fem::core