#pragma once

#include "analysis/l2_error.hpp"
#include "analysis/error_estimate.hpp"
#include "analysis/jump_indicator.hpp"
#include "analysis/refine.hpp"

#include "core/mesh.hpp"
#include "core/mesh_utils.hpp"

#include "problems/problem.hpp"
#include "problems/poisson_1d.hpp"

#include "assembly/assembler.hpp"
#include "assembly/sparsity.hpp"

#include "boundary/dirichlet_bc.hpp"
#include "boundary/neumann_bc.hpp"

#include "linalg/vector.hpp"
#include "linalg/solver.hpp"
#include "linalg/sparse_matrix.hpp"
#include "linalg/checks.hpp"

#include "solve/driver.hpp"

#include "io/output.hpp"
#include "io/meta.hpp"

#include "discretization/element/finite_element.hpp"
#include "discretization/element/lagrange_p1_1d.hpp"
#include "discretization/element/lagrange_p2_1d.hpp"
#include "discretization/element/lagrange_p1_2d.hpp"
#include "discretization/quadrature/quadrature_rule.hpp"
#include "discretization/quadrature/gauss_legendre_1d.hpp"