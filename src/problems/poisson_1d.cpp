#include "fem/problems/poisson_1d.hpp"

#include <stdexcept>
#include <utility>

namespace fem::problems {
/** 
Poisson1D::Poisson1D(ScalarFunction rhs)
    : m_rhs(std::move(rhs)),
      m_diffusion([](auto) { return 1.0; }),
      m_reaction([](auto) { return 0.0; }),
      m_exact([](auto) { return 0.0; }),
      m_name("Poisson1D"),
      m_meta({})
{
    if (!m_rhs) {
        throw std::invalid_argument("Poisson1D: rhs function must not be empty.");
    }
}

Poisson1D::Poisson1D(ScalarFunction rhs,
                     ScalarFunction diffusion)
    : m_rhs(std::move(rhs)),
      m_diffusion(std::move(diffusion)),
      m_reaction([](auto) { return 0.0; }),
      m_exact([](auto) { return 0.0; }),
      m_name("Poisson1D"),
      m_meta({})
{
    if (!m_rhs) {
        throw std::invalid_argument("Poisson1D: rhs function must not be empty.");
    }
    if (!m_diffusion) {
        throw std::invalid_argument("Poisson1D: diffusion function must not be empty.");
    }
}

Poisson1D::Poisson1D(ScalarFunction rhs,
                     ScalarFunction diffusion,
                     ScalarFunction reaction)
    : m_rhs(std::move(rhs)),
      m_diffusion(std::move(diffusion)),
      m_reaction(std::move(reaction)),
      m_exact([](auto) { return 0.0; }),
      m_name("Poisson1D"),
      m_meta({})
{
    if (!m_rhs) {
        throw std::invalid_argument("Poisson1D: rhs function must not be empty.");
    }
    if (!m_diffusion) {
        throw std::invalid_argument("Poisson1D: diffusion function must not be empty.");
    }
    if (!m_reaction) {
        throw std::invalid_argument("Poisson1D: reaction function must not be empty.");
    }
}
*/
Poisson1D::Poisson1D(ScalarFunction rhs,
                     ScalarFunction diffusion,
                     ScalarFunction reaction,
                     ScalarFunction exact,
                     std::string name,
                     fem::io::Meta meta)
    : m_rhs(std::move(rhs)),
      m_diffusion(std::move(diffusion)),
      m_reaction(std::move(reaction)),
      m_exact(std::move(exact)),
      m_name(std::move(name)),
      m_meta(std::move(meta))
{
    if (!m_rhs) {
        throw std::invalid_argument("Poisson1D: rhs function must not be empty.");
    }
    if (!m_diffusion) {
        throw std::invalid_argument("Poisson1D: diffusion function must not be empty.");
    }
    if (!m_reaction) {
        throw std::invalid_argument("Poisson1D: reaction function must not be empty.");
    }
    if (!m_exact) {
        throw std::invalid_argument("Poisson1D: exact function must not be empty.");
    }
}

double Poisson1D::rhs(double x) const {
    return m_rhs(x);
}

double Poisson1D::diffusion(double x) const {
    return m_diffusion(x);
}

double Poisson1D::reaction(double x) const {
    return m_reaction(x);
}

double Poisson1D::exact(double x) const {
    return m_exact(x);
}

std::string Poisson1D::name() const {
    return m_name;
}

fem::io::Meta Poisson1D::meta() const {
    return m_meta;
}

} // namespace fem::problems