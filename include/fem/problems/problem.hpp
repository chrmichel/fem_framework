#pragma once

#include <string>
#include <variant>
#include <map>
#include "fem/io/meta.hpp"

namespace fem {

namespace boundary {
}

namespace problems {
/**
 * Abstract PDE description for 1D elliptic problems:
 *
 *   - (a(x) u')' + c(x) u = f(x)
 *
 * Phase 2 design:
 *   Only rhs() is mandatory.
 *   All other coefficients have sensible defaults.
 */
class Problem {
public:
    virtual ~Problem() = default;

    // PDE-Koeffizienten (Phase 1: 1D Poisson)
    virtual double rhs(double x) const = 0;        // f(x)
    virtual double diffusion(double x) const {
        return 1.0; // a(x), default to constant diffusion
    }

    /// Reaction coefficient c(x)  (default: 0)
    virtual double reaction(double /*x*/) const {
        return 0.0;
    }

    /// Optional exact solution (default: 0)
    /// Override only in tests or manufactured solutions.
    virtual double exact(double /*x*/) const {
        return 0.0;
    }

    // NEW: minimal metadata hooks
    virtual std::string name() const {return "Problem";};           // e.g. "Poisson1D"
    virtual fem::io::Meta meta() const {return {};};         // problem-specific metadata (labels, etc.)
};

} // namespace problems
} // namespace fem