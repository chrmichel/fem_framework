# FEM Framework

A modular, dimension-generic C++20 Finite Element Method framework for elliptic PDEs.

---

## Overview

This framework solves second-order elliptic boundary value problems of the form

```
-(a(x) u')' + c(x) u = f(x)     in 1D
-div(a(x,y) grad u) + c(x,y) u = f(x,y)     in 2D
```

with Dirichlet and Neumann boundary conditions, using the finite element method.

The design is template-based on the spatial dimension `Dim`, allowing the same
assembly pipeline to work in 1D, 2D, and 3D without code duplication.

---

## Features

### Discretization
- `Mesh<Dim>` — dimension-generic mesh with explicit connectivity table
- `LagrangeP1_1D` — linear P1 element on reference segment `[-1, 1]`
- `LagrangeP2_1D` — quadratic P2 element on reference segment `[-1, 1]`
- `LagrangeP1_2D` — linear P1 triangle element on reference triangle `{(0,0),(1,0),(0,1)}`
- `GaussLegendre1D` — Gauss-Legendre quadrature on `[-1,1]` for n = 1..5
- `TriangleQuadrature` — Gauss quadrature on reference triangle for n = 1, 3

### Assembly
- `Assembler<Dim>` — generic assembler using `if constexpr` for dimension-specific Jacobi logic
- `build_sparse_matrix<Dim>` — builds CRS sparsity pattern from mesh and element
- 1D: scalar Jacobian, isoparametric mapping
- 2D: 2x2 Jacobian matrix, inverse via cofactor formula

### Boundary Conditions
- `BoundaryCondition<Dim>` — abstract base using `Point`-valued functions
- `DirichletBC<Dim>` — inhomogeneous Dirichlet via row elimination with RHS correction
  - 1D: applies to left and/or right boundary node
  - 2D: applies to all nodes in `boundary_node_ids()`
- `NeumannBC<Dim>` — natural boundary condition via RHS flux term
  - supports optional left and/or right boundary

### Linear Algebra
- `SparseMatrix` — Compressed Row Storage (CRS) with bounds-checked access
- `Vector` — dense vector with bounds checking
- `solve(SparseMatrix, Vector)` — Conjugate Gradient solver

### Solver
- `Driver<Dim>` — high-level facade: builds matrix, assembles, applies BCs, solves

### Analysis
- `jump_indicator` — elementwise jump indicator `h_e * |[a u']|` for adaptive refinement
- `refine` — local h-refinement via midpoint insertion, threshold `eta_e > theta * max(eta)`
- `l2_error<Dim>` — L2 error against exact solution, works for P1 and P2

### I/O
- `write_results<Dim>` — writes `solution.csv` and `meta.json` to a timestamped directory
- CSV header adapts automatically to dimension: `node_id, x [,y ,z], u`
- `meta.json` records problem name, element type, quadrature, boundary conditions, mesh size

---

## Architecture

```
apps/
  poisson1d/          — example application

include/fem/
  core/               — Mesh<Dim>, mesh_utils
  problems/           — Problem base class
  discretization/
    element/          — FiniteElement interface, LagrangeP1_1D, LagrangeP2_1D, LagrangeP1_2D
    quadrature/       — QuadratureRule interface, GaussLegendre1D, TriangleQuadrature
  assembly/           — Assembler<Dim>, sparsity pattern
  boundary/           — BoundaryCondition<Dim>, DirichletBC<Dim>, NeumannBC<Dim>
  linalg/             — SparseMatrix, Vector, solver
  solve/              — Driver<Dim>
  analysis/           — jump_indicator, refine, l2_error
  io/                 — write_results<Dim>, Meta

src/                  — implementations and explicit template instantiations
tests/                — test suite
```

### Design Principles

- Strict separation of concerns: mesh, PDE, discretization, assembly, and solver are independent layers
- No global state
- Dependency injection throughout — `Assembler` receives `Mesh`, `Problem`, `FiniteElement`, and `QuadratureRule` as references
- `BoundaryCondition` modifies the assembled system only, never the mesh or problem
- Template parameter `Dim` is the single source of truth for spatial dimension
- Explicit template instantiations (`template class Mesh<1>; template class Mesh<2>;`) keep compilation times manageable

---

## Build

### Requirements

- CMake ≥ 3.20
- C++20 compatible compiler (Clang 14+, GCC 12+, AppleClang 14+)

### Build and test

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

---

## Usage

### 1D Poisson with Dirichlet boundary conditions

```cpp
#include <fem/fem.hpp>
#include <cmath>

int main() {
    // Mesh: 80 uniform elements on [0, 1]
    std::vector<double> xs;
    for (std::size_t i = 0; i <= 80; ++i)
        xs.push_back(i / 80.0);
    fem::core::Mesh<1> mesh(
        [&]{ std::vector<fem::core::Mesh<1>::Point> pts; for (auto x : xs) pts.push_back({x}); return pts; }()
    );

    // Problem: -u'' = pi^2 sin(pi x), exact solution sin(pi x)
    fem::problems::Poisson1D problem(
        [](double x){ return M_PI * M_PI * std::sin(M_PI * x); }
    );

    // Boundary conditions: u(0) = 0, u(1) = 0
    fem::boundary::DirichletBC<1> bc(
        [](auto) { return 0.0; },
        [](auto) { return 0.0; }
    );

    auto fe   = fem::discretization::element::LagrangeP1_1D();
    auto quad = fem::discretization::quadrature::GaussLegendre1D(2);

    auto u = fem::Driver<1>::solve(mesh, problem, fe, quad, {&bc});

    fem::io::write_results(mesh, u, problem, fe, quad, {&bc});
}
```

### 1D mixed Dirichlet-Neumann

```cpp
// u(0) = 0 (Dirichlet), u'(1) = g (Neumann)
fem::boundary::DirichletBC<1> dirichlet(
    [](auto) { return 0.0; },
    std::nullopt
);
fem::boundary::NeumannBC<1> neumann(
    std::nullopt,
    [](auto) { return 2.0; }
);

auto u = fem::Driver<1>::solve(mesh, problem, fe, quad, {&neumann, &dirichlet});
```

### 2D Poisson on a triangle mesh

```cpp
// 16x16 triangle mesh on [0,1]^2
auto mesh = fem::core::create_triangle_mesh_2d(0.0, 1.0, 0.0, 1.0, 16, 16);

// Problem: -Delta u = 2 pi^2 sin(pi x) sin(pi y)
struct MyProblem : fem::problems::Problem {
    double rhs(double x, double y) const override {
        return 2.0 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
    }
};

MyProblem problem;
fem::boundary::DirichletBC<2> bc(
    [](auto) { return 0.0; },
    std::nullopt
);

auto fe   = fem::discretization::element::LagrangeP1_2D();
auto quad = fem::discretization::quadrature::TriangleQuadrature(3);

auto u = fem::Driver<2>::solve(mesh, problem, fe, quad, {&bc});
```

### Adaptive refinement (1D)

```cpp
auto est  = fem::analysis::jump_indicator(mesh, u, problem);
auto mesh2 = fem::analysis::refine(mesh, est, 0.5);  // refine elements with eta > 0.5 * max(eta)
auto u2   = fem::Driver<1>::solve(mesh2, problem, fe, quad, {&bc});
```

---

## Test Suite

| Test | What it verifies |
|---|---|
| `mesh_test` | Mesh construction, connectivity, boundary nodes |
| `mesh2d_test` | 2D triangle mesh, node coordinates, connectivity |
| `linear_dirichlet_test` | P1 exactness for linear solutions |
| `convergence_test` | P1 L2 convergence rate > 1.5 |
| `l2_error_test` | L2 error decreases with refinement |
| `reaction_test` | Reaction term changes the solution |
| `neumann_bc_test` | Mixed Dirichlet-Neumann, L2 convergence |
| `p2_element_test` | P2 shape functions, partition of unity, derivatives |
| `p2_convergence_test` | P2 L2 convergence rate > 1.8 |
| `output_test` | CSV and JSON output format |
| `sparse_matrix_test` | CRS access, set_zero, out-of-pattern throws |
| `sparsity_test` | Pattern correctness for P1 and P2 |
| `adaptive_refinement_test` | Jump indicator, refinement, error decrease |
| `p1_2d_element_test` | 2D shape functions, partition of unity, derivatives |
| `triangle_quadrature_test` | Weights sum, polynomial exactness |
| `poisson2d_test` | 2D assembler: positive diagonal, zero row sums |
| `poisson2d_convergence_test` | 2D P1 L2 convergence rate > 1.5 |

---

## Author

Christian Michel

## License

MIT License