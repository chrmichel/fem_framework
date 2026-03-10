#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace fem::test {

inline void fail(const std::string& msg) {
    std::cerr << "[TEST FAIL] " << msg << "\n";
    std::exit(1);
}

inline void require(bool cond, const std::string& msg) {
    if (!cond) fail(msg);
}

inline void require_close(double a, double b, double tol, const std::string& msg) {
    if (std::fabs(a - b) > tol) {
        std::cerr << "[TEST FAIL] " << msg << " | a=" << a << " b=" << b << " tol=" << tol << "\n";
        std::exit(1);
    }
}

} // namespace fem::test