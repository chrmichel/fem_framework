#include "fem/core/mesh.hpp"

#include <stdexcept>

namespace fem {

Mesh::Mesh(double a, double b, Index elements)
    : dim_(1), a_(a), b_(b), elements_(elements)
{
    if (elements_ == 0) {
        throw std::invalid_argument("Mesh: elements must be > 0");
    }
    if (!(b_ > a_)) {
        throw std::invalid_argument("Mesh: require b > a");
    }

    const double h = (b_ - a_) / static_cast<double>(elements_);

    // Nodes
    nodes_.resize(elements_ + 1);
    for (Index i = 0; i < nodes_.size(); ++i) {
        const double x = a_ + static_cast<double>(i) * h;
        nodes_[i] = Point{x, 0.0, 0.0};
    }

    // Connectivity
    elements_conn_.resize(elements_);
    for (Index e = 0; e < elements_; ++e) {
        elements_conn_[e] = {e, e + 1};
    }
}

Mesh::Index Mesh::dimension() const noexcept { return dim_; }
Mesh::Index Mesh::n_nodes() const noexcept { return nodes_.size(); }
Mesh::Index Mesh::n_elements() const noexcept { return elements_conn_.size(); }

const Mesh::Point& Mesh::point(Index i) const
{
    if (i >= nodes_.size()) {
        throw std::out_of_range("Mesh::point: node index out of range");
    }
    return nodes_[i];
}

double Mesh::node(Index i) const
{
    return point(i).x; // 1D convenience
}

std::span<const Mesh::Index> Mesh::cell_nodes(Index e) const
{
    if (e >= elements_conn_.size()) {
        throw std::out_of_range("Mesh::cell_nodes: element index out of range");
    }
    // span über die 2 Node-IDs des 1D-Elements
    return std::span<const Index>(elements_conn_[e].data(), elements_conn_[e].size());
}

std::array<Mesh::Index, 2> Mesh::element(Index e) const
{
    if (e >= elements_conn_.size()) {
        throw std::out_of_range("Mesh::element: element index out of range");
    }
    return elements_conn_[e];
}

std::vector<Mesh::Index> Mesh::boundary_nodes() const
{
    // Phase 1: Intervall hat genau zwei Randknoten
    return { left_boundary_node(), right_boundary_node() };
}

Mesh::Index Mesh::left_boundary_node() const noexcept { return 0; }
Mesh::Index Mesh::right_boundary_node() const noexcept { return n_nodes() - 1; }

} // namespace fem