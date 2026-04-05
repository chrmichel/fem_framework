#pragma once

#include "boundary_condition.hpp"

#include "fem/core/mesh.hpp"
#include "fem/linalg/sparse_matrix.hpp"
#include "fem/linalg/vector.hpp"
#include "fem/linalg/checks.hpp"
#include <functional>
#include <cstddef>
#include <optional>

namespace fem {
class Mesh;

namespace boundary {

template<int Dim>
class DirichletBC final : public BoundaryCondition<Dim> {
public:
    using Function = std::function<double(double)>;

    DirichletBC(std::optional<Function> g_left, std::optional<Function> g_right)
        : g_left_(std::move(g_left)), g_right_(std::move(g_right))
    {
        if (!g_left_ && !g_right_) {
            throw std::invalid_argument("DirichletBC: at least one boundary function must be set.");
        }
    }

    void apply(linalg::SparseMatrix& A,
               linalg::Vector& b,
                const core::Mesh<Dim>& mesh) const override{
        linalg::require_square(A);
        linalg::require_same_size(A, b);

        if (g_left_) {
            const std::size_t iL = mesh.left_boundary_node();
            enforce_node(iL, (*g_left_)(mesh.node(iL)[0]), A, b);
        }

        if (g_right_) {
            const std::size_t iR = mesh.right_boundary_node();
            enforce_node(iR, (*g_right_)(mesh.node(iR)[0]), A, b);
        }
    }
    
    virtual std::string name() const override {
        return "DirichletBC";
    }

private:
    std::optional<Function> g_left_;
    std::optional<Function> g_right_;

    static void enforce_node(std::size_t node,
                             double value,
                             linalg::SparseMatrix& A,
                             linalg::Vector& b) {
        const std::size_t n = A.rows();

        for (std::size_t j = 0; j < n; ++j) {
            // nur wenn (j, node) im Pattern existiert
            try {
                b(j) -= A(j, node) * value;
                A(j, node) = 0.0;
            } catch (const std::out_of_range&) {
                // Eintrag nicht im Pattern — nichts zu tun
            }
        }

        // Zeile node nullen
        for (std::size_t k = A.row_ptr()[node]; k < A.row_ptr()[node + 1]; ++k) {
            A.values_mutable()[k] = 0.0;
        }

        A(node, node) = 1.0;
        b(node) = value;
    }
};

} // namespace boundary
} // namespace fem