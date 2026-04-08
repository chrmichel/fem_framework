#pragma once
#include "fem/core/mesh.hpp"
#include "fem/core/mesh_impl.hpp"
#include "fem/linalg/sparse_matrix.hpp"
#include "fem/discretization/element/finite_element.hpp"

#include <set>
#include <vector>

namespace fem::assembly {

template<int Dim>
linalg::SparseMatrix build_sparse_matrix(
    const core::Mesh<Dim>& mesh,
    const discretization::element::FiniteElement& fe)
{
    const std::size_t n     = mesh.n_nodes();
    const std::size_t ndofs = fe.num_dofs();

    std::vector<std::set<std::size_t>> pattern(n);

    for (std::size_t e = 0; e < mesh.n_elements(); ++e) {
        for (std::size_t i = 0; i < ndofs; ++i) {
            const std::size_t gi = mesh.element_node(e, i);
            for (std::size_t j = 0; j < ndofs; ++j) {
                const std::size_t gj = mesh.element_node(e, j);
                pattern[gi].insert(gj);
            }
        }
    }

    std::vector<std::size_t> row_ptr(n + 1, 0);
    std::vector<std::size_t> col_idx;

    for (std::size_t i = 0; i < n; ++i) {
        row_ptr[i + 1] = row_ptr[i] + pattern[i].size();
        for (std::size_t j : pattern[i])
            col_idx.push_back(j);
    }

    return linalg::SparseMatrix(n, std::move(row_ptr), std::move(col_idx));
}

} // namespace fem::assembly