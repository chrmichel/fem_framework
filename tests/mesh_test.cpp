#include "tests/test_assert.hpp"
#include "fem/core/mesh.hpp"

#include <vector>

int main() {
    fem::core::Mesh mesh(std::vector<double>{0.0, 0.5, 1.0});

    fem::test::require(mesh.n_nodes() == 3, "mesh.n_nodes() == 3");
    fem::test::require(mesh.n_elements() == 2, "mesh.n_elements() == 2");

    fem::test::require(mesh.node(0) == 0.0, "mesh.node(0) == 0");
    fem::test::require(mesh.node(1) == 0.5, "mesh.node(1) == 0.5");
    fem::test::require(mesh.node(2) == 1.0, "mesh.node(2) == 1");

    fem::test::require(mesh.element_node(0,0) == 0, "element_node(0,0) == 0");
    fem::test::require(mesh.element_node(0,1) == 1, "element_node(0,1) == 1");
    fem::test::require(mesh.element_node(1,0) == 1, "element_node(1,0) == 1");
    fem::test::require(mesh.element_node(1,1) == 2, "element_node(1,1) == 2");

    // optional boundary helpers, only enable if they exist
    // fem::test::require(mesh.left_boundary_node() == 0, "left_boundary_node == 0");
    // fem::test::require(mesh.right_boundary_node() == mesh.n_nodes()-1, "right_boundary_node == last");

    return 0;
}