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
    fem::test::require(mesh.left_boundary_node() == 0, "left_boundary_node == 0");
    fem::test::require(mesh.right_boundary_node() == mesh.n_nodes()-1, "right_boundary_node == last");

    // P1-Mesh mit expliziter Konnektivität — muss identisch zu auto-generiertem verhalten
    std::vector<std::vector<std::size_t>> conn = {{0,1},{1,2},{2,3}};
    fem::core::Mesh mesh2(std::vector<double>{0.0, 0.25, 0.75, 1.0}, conn);

    fem::test::require(mesh2.n_elements() == 3, "explicit conn: n_elements == 3");
    fem::test::require(mesh2.element_node(0,0) == 0, "explicit conn: element_node(0,0)");
    fem::test::require(mesh2.element_node(2,1) == 3, "explicit conn: element_node(2,1)");

    // Index out of range soll werfen
    bool threw = false;
    try {
        fem::core::Mesh bad(std::vector<double>{0.0, 1.0}, {{0, 5}});
    } catch (const std::exception&) {
        threw = true;
    }
    fem::test::require(threw, "explicit conn: out-of-range index throws");

    return 0;
}