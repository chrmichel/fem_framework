#pragma once

#include <array>
#include <vector>
#include <cstddef>

namespace fem::core {

template<int Dim>
class Mesh {
public:
    using Point = std::array<double, Dim>;

    explicit Mesh(std::vector<Point> nodes);

    Mesh(std::vector<Point> nodes,
         std::vector<std::vector<std::size_t>> connectivity);

    std::size_t n_nodes()    const;
    std::size_t n_elements() const;
    std::size_t nodes_per_element(std::size_t e) const;

    Point node(std::size_t i) const;

    std::size_t element_node(std::size_t e,
                             std::size_t local_index) const;

    std::size_t left_boundary_node()  const;
    std::size_t right_boundary_node() const;
    std::array<std::size_t, 2> boundary_nodes() const;

    void set_boundary_nodes(std::size_t left, std::size_t right);

private:
    std::vector<Point> m_nodes;
    std::vector<std::vector<std::size_t>> m_connectivity;
    std::size_t m_left_boundary_node{0};
    std::size_t m_right_boundary_node{0};

    void validate_nodes()       const;
    void validate_connectivity() const;
};

extern template class Mesh<1>;
extern template class Mesh<2>;
extern template class Mesh<3>;

using Mesh1D = Mesh<1>;
using Mesh2D = Mesh<2>;
using Mesh3D = Mesh<3>;

} // namespace fem::core