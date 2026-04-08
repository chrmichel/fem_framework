#pragma once

#include <string>
#include <variant>
#include <map>
#include "fem/io/meta.hpp"

namespace fem {

namespace problems {

class Problem {
public:
    virtual ~Problem() = default;

    // --- 1D ---
    virtual double rhs(double x) const { return 0.0; }
    virtual double diffusion(double x) const { return 1.0; }
    virtual double reaction(double x) const { return 0.0; }
    virtual double exact(double x) const { return 0.0; }

    // --- 2D --- defaults delegieren an 1D (für Kompatibilität)
    virtual double rhs(double x, double y) const {
        return rhs(x);
    }
    virtual double diffusion(double x, double y) const {
        return diffusion(x);
    }
    virtual double reaction(double x, double y) const {
        return reaction(x);
    }
    virtual double exact(double x, double y) const {
        return exact(x);
    }

    virtual std::string name() const { return "Problem"; }
    virtual fem::io::Meta meta() const { return {}; }
};

} // namespace problems
} // namespace fem