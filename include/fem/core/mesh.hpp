#pragma once
#include <cstddef>
#include <vector>
#include <array>
#include <span>

namespace fem
{
class Mesh{
public:
    using Index = std::size_t;

    //Points
    struct Point{
        double x{0.0}, y{0.0}, z{0.0};
        double operator[](Index i) const noexcept {
            return (i==0)?x:(i==1)?y:z;
        }
    };
    //ctor
    Mesh(double a, double b, Index elements);
    //meta
    Index dimension() const noexcept;
    Index n_nodes() const noexcept;
    Index n_elements() const noexcept;
    //geometry
    const Point& point(Index i) const;
    double node(Index i) const; //1D shorthand
    //topology
    std::span<const Index> cell_nodes(Index e) const;
    std::array<Index,2> element(Index e) const;
    //boundary
    std::vector<Index> boundary_nodes() const;
    Index left_boundary_node() const noexcept;
    Index right_boundary_node() const noexcept;

private:
    Index dim_{1};
    double a_{0.0};
    double b_{0.0};
    Index elements_{0};

    std::vector<Point> nodes_;
    //connectivity
    std::vector<std::array<Index,2>> elements_conn_;
};
} // namespace fem