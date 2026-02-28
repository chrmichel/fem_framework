#pragma once

#include "fem/problems/problem.hpp"
#include "fem/linalg/matrix.hpp"
#include "fem/linalg/vector.hpp"

namespace fem::assembly {

class Assembler {
public:
    // Baut A und b (ohne BCs)
    fem::linalg::Matrix assemble_stiffness(const fem::problems::Problem& problem) const;
    fem::linalg::Vector assemble_load(const fem::problems::Problem& problem) const;

    // Optional convenience: baut A,b und wendet alle BCs an
    void assemble_system(const fem::problems::Problem& problem,
                         fem::linalg::Matrix& A,
                         fem::linalg::Vector& b) const;

private:
    // P1-Shape-Funktionen auf Referenzelement [-1,1]
    static double phi(int a, double xi);
    static double dphi_dxi(int a);
};

} // namespace fem::assembly