#pragma once
#include <vector>

namespace fem::analysis {

struct ErrorEstimate {
    std::vector<double> indicators;  // ein Wert pro Element
    double global_error;             // sqrt(sum eta_e^2)
};

} // namespace fem::analysis