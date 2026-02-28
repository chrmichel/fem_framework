#include <fem/fem.hpp>
#include <iostream>

int main()
{
    fem::problems::Poisson1D problem(
        0.0, 1.0, 80,
        [](double x){ return 1.0; },   // rhs
        [](double x){ return 1.0; },   // diffusion
        [](double){ return 2.0; },     // g_left
        [](double){ return -1.0; }     // g_right
    );

    fem::solve::Driver driver;
    auto u = driver.solve(problem);

    std::cout << "u size = " << u.size() << "\n";
}