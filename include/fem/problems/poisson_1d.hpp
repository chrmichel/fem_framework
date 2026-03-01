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
              Function g_right,
              std::string rhs_label = "custom",
              std::string diffusion_label = "custom",
              std::string g_left_label = "custom",
              std::string g_right_label = "custom")
        : mesh_(std::make_shared<fem::Mesh>(a, b, elements)),
          rhs_(std::move(rhs)),
          diffusion_(std::move(diffusion)),
          rhs_label_(std::move(rhs_label)),
          diffusion_label_(std::move(diffusion_label)),
          g_left_label_(std::move(g_left_label)),
          g_right_label_(std::move(g_right_label))
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
    std::string name() const override { return "Poisson1D"; }

    fem::io::Meta meta() const override {
        fem::io::Meta m;
        m["rhs"] = rhs_label_;
        m["diffusion"] = diffusion_label_;
        m["dirichlet_left"] = g_left_label_;
        m["dirichlet_right"] = g_right_label_;
        m["bc_type"] = std::string("Dirichlet");
        return m;
    }

private:
    std::shared_ptr<fem::Mesh> mesh_;
    Function rhs_;
    Function diffusion_;
    std::string rhs_label_;
    std::string diffusion_label_;
    std::string g_left_label_;
    std::string g_right_label_;
    std::vector<std::shared_ptr<const fem::boundary::BoundaryCondition>> bcs_;
};

} // namespace fem::problems