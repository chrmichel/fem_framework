#pragma once

#include "analysis/l2_error.hpp"

#include "core/mesh.hpp"

#include "problems/problem.hpp"
#include "problems/poisson_1d.hpp"

#include "assembly/assembler.hpp"

#include "boundary/dirichlet_bc.hpp"

#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"
#include "linalg/solver.hpp"

#include "solve/driver.hpp"

#include "io/output.hpp"
#include "io/meta.hpp"

#include "discretization/element/finite_element.hpp"
#include "discretization/element/lagrange_p1_1d.hpp"
#include "discretization/quadrature/quadrature_rule.hpp"
#include "discretization/quadrature/gauss_legendre_1d.hpp"