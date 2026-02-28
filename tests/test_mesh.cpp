#include <cassert>
#include <iostream>
#include <vector>
#include <fem/core/mesh.hpp>

int main()
{
    fem::Mesh mesh(0.0, 1.0, 10);

    assert(mesh.dimension() == 1);
    assert(mesh.n_elements() == 10);
    assert(mesh.n_nodes() == 11);

    // boundary nodes
    auto bn = mesh.boundary_nodes();
    assert(bn.size() == 2);
    assert(bn[0] == 0);
    assert(bn[1] == mesh.n_nodes() - 1);

    assert(mesh.left_boundary_node() == 0);
    assert(mesh.right_boundary_node() == mesh.n_nodes() - 1);

    // connectivity check
    for (std::size_t e = 0; e < mesh.n_elements(); ++e) {
        auto el = mesh.element(e);
        assert(el[0] == e);
        assert(el[1] == e + 1);

        // cell_nodes span should match
        auto cn = mesh.cell_nodes(e);
        assert(cn.size() == 2);
        assert(cn[0] == e);
        assert(cn[1] == e + 1);
    }

    std::cout << "MeshTest passed.\n";
    return 0;
}