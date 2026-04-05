#include "tests/test_assert.hpp"
#include "fem/core/mesh_utils.hpp"
#include "fem/core/mesh_impl.hpp"

int main() {
    // 2x2 Gitter -> 4 Rechtecke -> 8 Dreiecke, 9 Knoten
    const auto mesh = fem::core::create_triangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 2, 2);

    fem::test::require(mesh.n_nodes()    == 9, "n_nodes == 9");
    fem::test::require(mesh.n_elements() == 8, "n_elements == 8");
    fem::test::require(mesh.nodes_per_element(0) == 3, "nodes_per_element == 3");

    // Erster Knoten ist (0,0)
    fem::test::require_close(mesh.node(0)[0], 0.0, 1e-14, "node(0).x == 0");
    fem::test::require_close(mesh.node(0)[1], 0.0, 1e-14, "node(0).y == 0");

    // Letzter Knoten ist (1,1)
    fem::test::require_close(mesh.node(8)[0], 1.0, 1e-14, "node(8).x == 1");
    fem::test::require_close(mesh.node(8)[1], 1.0, 1e-14, "node(8).y == 1");

    // Konnektivität: erstes Dreieck hat 3 gültige Indizes
    for (std::size_t local = 0; local < 3; ++local) {
        const std::size_t idx = mesh.element_node(0, local);
        fem::test::require(idx < mesh.n_nodes(),
            "element_node(0," + std::to_string(local) + ") in range");
    }

    // Alle Elemente haben 3 Knoten
    for (std::size_t e = 0; e < mesh.n_elements(); ++e) {
        fem::test::require(mesh.nodes_per_element(e) == 3,
            "element " + std::to_string(e) + " has 3 nodes");
    }

    // 1x1 Gitter -> 2 Dreiecke, 4 Knoten
    const auto mesh2 = fem::core::create_triangle_mesh_2d(
        0.0, 2.0, 0.0, 3.0, 1, 1);

    fem::test::require(mesh2.n_nodes()    == 4, "1x1: n_nodes == 4");
    fem::test::require(mesh2.n_elements() == 2, "1x1: n_elements == 2");

    // Knotenkoordinaten prüfen
    fem::test::require_close(mesh2.node(1)[0], 2.0, 1e-14, "1x1: node(1).x == 2");
    fem::test::require_close(mesh2.node(2)[1], 3.0, 1e-14, "1x1: node(2).y == 3");

    return 0;
}