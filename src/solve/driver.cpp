#include "fem/solve/driver.hpp"

namespace fem::solve {

fem::linalg::Vector Driver::solve(const fem::problems::Problem& problem) const
{
    fem::linalg::Matrix A;
    fem::linalg::Vector b;

    // Assembler baut System und wendet BCs an
    assembler_.assemble_system(problem, A, b);

    return fem::linalg::solve(A, b);
}

Result Driver::solve_with_system(const fem::problems::Problem& problem) const
{
    Result r;

    assembler_.assemble_system(problem, r.A, r.b);
    r.u = fem::linalg::solve(r.A, r.b);

    return r;
}

} // namespace fem::solve