#pragma once

#include <memory>
#include <vector>

namespace fem {

class Mesh;

namespace boundary {
class BoundaryCondition;
}

namespace problems {

class Problem {
public:
    virtual ~Problem() = default;

    // Geometrie / Diskretisierungsraum
    virtual const fem::Mesh& mesh() const = 0;

    // PDE-Koeffizienten (Phase 1: 1D Poisson)
    virtual double rhs(double x) const = 0;        // f(x)
    virtual double diffusion(double x) const = 0;  // a(x)

    // Randbedingungen als Strategie-Objekte
    virtual std::vector<std::shared_ptr<const fem::boundary::BoundaryCondition>>
    boundary_conditions() const = 0;
};

} // namespace problems
} // namespace fem