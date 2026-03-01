#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <variant>

namespace fem::io {

using MetaValue = std::variant<std::string, double, std::int64_t, bool>;
using Meta      = std::map<std::string, MetaValue>;

} // namespace fem::io