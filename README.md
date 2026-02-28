# FEMFramework

A modular C++20 Finite Element Method (FEM) framework.

Phase 1 completed:  
1D Poisson equation with inhomogeneous Dirichlet boundary conditions, full test coverage, clean architecture, and extensible design.

---

## Features (Phase 1)

- 1D uniform mesh (2D-ready internal structure)
- Linear P1 finite elements
- Variable diffusion coefficient \( a(x) \)
- Arbitrary right-hand side \( f(x) \)
- Inhomogeneous Dirichlet boundary conditions
- Polymorphic `Problem` interface
- BoundaryCondition strategy pattern
- Dense Matrix/Vector with internal safety checks
- Gaussian elimination solver with partial pivoting
- High-level `solve::Driver` facade
- Complete test suite including convergence and L2 error validation

---

## Mathematical Model

We solve the 1D Poisson problem

\[
-(a(x) u')' = f(x), \quad x \in (a,b)
\]

with Dirichlet boundary conditions

\[
u(a) = g_L, \quad u(b) = g_R.
\]

Discretization:
- Linear (P1) finite elements
- 2-point Gauss quadrature
- Consistent assembly
- Proper inhomogeneous Dirichlet elimination

---
### Design Principles

- Strict separation of concerns
- No global state
- Mesh independent of PDE
- Assembly independent of boundary conditions
- Boundary conditions modify algebra only
- Solver independent of discretization
- Fully test-driven validation

---

## Build Instructions

### Requirements

- CMake ≥ 3.20
- C++20 compatible compiler (Clang / GCC / AppleClang)

### Build

```bash
cmake -S . -B build
cmake --build build
Run Tests
ctest --test-dir build --output-on-failure
```

## Example usage
```cpp
#include <fem/fem.hpp>

int main()
{
    fem::problems::Poisson1D problem(
        0.0, 1.0, 80,
        [](double x){ return 1.0; },   // rhs
        [](double x){ return 1.0; },   // diffusion
        [](double){ return 2.0; },     // left Dirichlet
        [](double){ return -1.0; }     // right Dirichlet
    );

    fem::solve::Driver driver;
    auto u = driver.solve(problem);
}
```

## License
MIT License

## Author
Christian Michel