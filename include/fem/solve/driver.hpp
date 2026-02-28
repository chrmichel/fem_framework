#pragma once

#include "fem/problems/problem.hpp"
#include "fem/assembly/assembler.hpp"
#include "fem/linalg/solver.hpp"
#include "fem/linalg/matrix.hpp"
#include "fem/linalg/vector.hpp"

namespace fem::solve {

// Optional: Ergebniscontainer, praktisch für Debug/Tests
struct Result {
    fem::linalg::Vector u;
    fem::linalg::Matrix A;
    fem::linalg::Vector b;
};

class Driver {
public:
    Driver() = default;

    // "High-level" solve: returns only u
    fem::linalg::Vector solve(const fem::problems::Problem& problem) const;

    // Debug-friendly: returns u and also A,b
    Result solve_with_system(const fem::problems::Problem& problem) const;

private:
    fem::assembly::Assembler assembler_;
};

} // namespace fem::solve