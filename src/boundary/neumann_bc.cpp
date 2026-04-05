#include "fem/boundary/neumann_bc.hpp"

#include <stdexcept>
#include <utility>

namespace fem::boundary {

template class NeumannBC<1>;
template class NeumannBC<2>;
template class NeumannBC<3>;

} // namespace fem::boundary