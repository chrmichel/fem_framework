#include "fem/assembly/assembler.hpp"

namespace fem::assembly {
    template class Assembler<1>;
    template class Assembler<2>;
    template class Assembler<3>;
}