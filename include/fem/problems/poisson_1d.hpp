#pragma once

#include "fem/problems/problem.hpp"
#include "fem/io/meta.hpp"

#include <functional>
#include <string>

namespace fem::problems {

/**
 * Generic 1D Poisson/reaction problem:
 *
 *   - (a(x) u')' + c(x) u = f(x)
 *
 * with user-provided coefficient functions.
 *
 * This class is intentionally lightweight and mesh-independent.
 */
class Poisson1D final : public Problem {
public:
    using ScalarFunction = std::function<double(double)>;

    // explicit Poisson1D(ScalarFunction rhs);

    // Poisson1D(ScalarFunction rhs,
    //           ScalarFunction diffusion);

    // Poisson1D(ScalarFunction rhs,
    //           ScalarFunction diffusion,
    //           ScalarFunction reaction);

    Poisson1D(ScalarFunction rhs,
                ScalarFunction diffusion = [](double) { return 1.0; },
                ScalarFunction reaction = [](double) { return 0.0; },
                ScalarFunction exact = [](double) { return 0.0; },
                std::string name = "Poisson1D",
                fem::io::Meta meta = {});

    double rhs(double x) const override;
    double diffusion(double x) const override;
    double reaction(double x) const override;
    double exact(double x) const override;

    std::string name() const override;
    fem::io::Meta meta() const override;

private:
    ScalarFunction m_rhs;
    ScalarFunction m_diffusion;
    ScalarFunction m_reaction;
    ScalarFunction m_exact;

    std::string m_name;
    fem::io::Meta m_meta;
};

} // namespace fem::problems