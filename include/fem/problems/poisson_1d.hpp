#pragma once

#include "problem.hpp"
#include "fem/core/mesh.hpp"
#include "fem/boundary/dirichlet_bc.hpp"

#include <functional>
#include <memory>
#include <utility>

namespace fem::problems {

class Poisson1D final : public Problem {
public:
    using Function = std::function<double(double)>;

    // Mesh wird intern gehalten (vollständiges Problemobjekt)
    Poisson1D(double a,
              double b,
              std::size_t elements,
              Function rhs,
              Function diffusion,
              Function g_left,
              Function g_right)
        : mesh_(std::make_shared<fem::Mesh>(a, b, elements)),
          rhs_(std::move(rhs)),
          diffusion_(std::move(diffusion))
    {
        // DirichletBC bekommt zwei Funktionen (links/rechts),
        // nutzt mesh.boundary_nodes() / left/right node intern.
        bcs_.push_back(std::make_shared<fem::boundary::DirichletBC>(
            std::move(g_left),
            std::move(g_right)
        ));
    }

    const fem::Mesh& mesh() const override
    {
        return *mesh_;
    }

    double rhs(double x) const override
    {
        return rhs_(x);
    }

    double diffusion(double x) const override
    {
        return diffusion_(x);
    }

    std::vector<std::shared_ptr<const fem::boundary::BoundaryCondition>>
    boundary_conditions() const override
    {
        // Upcast to interface type
        std::vector<std::shared_ptr<const fem::boundary::BoundaryCondition>> out;
        out.reserve(bcs_.size());
        for (const auto& bc : bcs_) out.push_back(bc);
        return out;
    }

private:
    std::shared_ptr<fem::Mesh> mesh_;
    Function rhs_;
    Function diffusion_;

    std::vector<std::shared_ptr<const fem::boundary::BoundaryCondition>> bcs_;
};

} // namespace fem::problems