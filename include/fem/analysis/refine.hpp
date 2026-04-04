#pragma once
#include "error_estimate.hpp"
#include "fem/core/mesh.hpp"

namespace fem::analysis {

/// Verfeinert alle Elemente mit eta_e > theta * max(eta)
/// theta in (0, 1], Standard: 0.5
core::Mesh1D refine(const core::Mesh1D& mesh,
                  const ErrorEstimate& estimate,
                  double theta = 0.5);

} // namespace fem::analysis