#include "fem/boundary/dirichlet_bc.hpp"

#include <stdexcept>
#include <utility>

namespace fem::boundary {

template class DirichletBC<1>;
template class DirichletBC<2>;
template class DirichletBC<3>;

} // namespace fem::boundary