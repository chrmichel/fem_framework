#pragma once
#include "error_estimate.hpp"
#include "fem/core/mesh.hpp"
#include "fem/linalg/vector.hpp"
#include "fem/problems/problem.hpp"

namespace fem::analysis {

ErrorEstimate jump_indicator(
    const core::Mesh& mesh,
    const linalg::Vector& u,
    const problems::Problem& problem);

} // namespace fem::analysis