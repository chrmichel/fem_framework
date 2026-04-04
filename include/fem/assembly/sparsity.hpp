#pragma once
#include "fem/core/mesh.hpp"
#include "fem/linalg/sparse_matrix.hpp"
#include "fem/discretization/element/finite_element.hpp"

namespace fem::assembly {

// Baut SparseMatrix mit korrektem Pattern aus Mesh + Element
linalg::SparseMatrix build_sparse_matrix(
    const core::Mesh& mesh,
    const discretization::element::FiniteElement& fe);

} // namespace fem::assembly