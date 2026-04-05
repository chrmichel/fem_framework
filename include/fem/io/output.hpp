#pragma once

#include "fem/core/mesh.hpp"
#include "fem/core/mesh_impl.hpp"
#include "fem/problems/problem.hpp"
#include "fem/io/meta.hpp"
#include "fem/discretization/element/finite_element.hpp"
#include "fem/discretization/quadrature/quadrature_rule.hpp"
#include "fem/boundary/boundary_condition.hpp"
#include "fem/linalg/vector.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

namespace fem::io {

namespace fs = std::filesystem;

namespace detail {

inline std::string make_timestamp()
{
    const auto now = std::chrono::system_clock::now();
    const auto t   = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    return oss.str();
}

inline std::string escape_json(const std::string& s)
{
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
        case '\"': out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\n': out += "\\n";  break;
        case '\r': out += "\\r";  break;
        case '\t': out += "\\t";  break;
        default:   out += c;      break;
        }
    }
    return out;
}

inline void write_meta_value(std::ostream& out, const MetaValue& value)
{
    std::visit(
        [&](const auto& v) {
            using T = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<T, std::string>) {
                out << "\"" << escape_json(v) << "\"";
            } else if constexpr (std::is_same_v<T, bool>) {
                out << (v ? "true" : "false");
            } else {
                out << v;
            }
        },
        value
    );
}

inline void write_meta_object(std::ostream& out,
                               const Meta& meta,
                               int indent = 2)
{
    out << "{\n";
    std::size_t count = 0;
    for (const auto& [key, value] : meta) {
        out << std::string(indent, ' ')
            << "\"" << escape_json(key) << "\": ";
        write_meta_value(out, value);
        ++count;
        if (count < meta.size()) out << ",";
        out << "\n";
    }
    out << std::string(indent - 2, ' ') << "}";
}

template<int Dim>
void write_solution_csv(const fs::path& file,
                        const core::Mesh<Dim>& mesh,
                        const linalg::Vector& solution)
{
    if (solution.size() != mesh.n_nodes()) {
        throw std::runtime_error(
            "write_solution_csv: solution size does not match mesh.n_nodes()");
    }

    std::ofstream out(file);
    if (!out) {
        throw std::runtime_error(
            "write_solution_csv: could not open file " + file.string());
    }

    // Header: node_id, x (,y ,z je nach Dim), u
    out << "node_id";
    if constexpr (Dim >= 1) out << ",x";
    if constexpr (Dim >= 2) out << ",y";
    if constexpr (Dim >= 3) out << ",z";
    out << ",u\n";

    for (std::size_t i = 0; i < mesh.n_nodes(); ++i) {
        out << i;
        for (int d = 0; d < Dim; ++d)
            out << "," << std::setprecision(16) << mesh.node(i)[d];
        out << "," << std::setprecision(16) << solution(i) << "\n";
    }
}

template<int Dim>
void write_meta_json(const fs::path& file,
                     const core::Mesh<Dim>& mesh,
                     const problems::Problem& problem,
                     const discretization::element::FiniteElement& fe,
                     const discretization::quadrature::QuadratureRule& quad,
                     const std::vector<const boundary::BoundaryCondition<Dim>*>& bcs)
{
    std::ofstream out(file);
    if (!out) {
        throw std::runtime_error(
            "write_meta_json: could not open file " + file.string());
    }

    std::string bc_names;
    for (const auto* bc : bcs) {
        if (!bc_names.empty()) bc_names += ", ";
        bc_names += bc->name();
    }

    out << "{\n";
    out << "  \"dim\": " << Dim << ",\n";
    out << "  \"problem\": \"" << escape_json(problem.name()) << "\",\n";
    out << "  \"problem_meta\": ";
    write_meta_object(out, problem.meta(), 4);
    out << ",\n";
    out << "  \"finite_element\": \"" << escape_json(fe.name()) << "\",\n";
    out << "  \"quadrature\": \"" << escape_json(quad.name()) << "\",\n";
    out << "  \"boundary_condition\": \"" << escape_json(bc_names) << "\",\n";
    out << "  \"n_nodes\": " << mesh.n_nodes() << ",\n";
    out << "  \"n_elements\": " << mesh.n_elements() << "\n";
    out << "}\n";
}

} // namespace detail

template<int Dim>
fs::path write_results(
    const core::Mesh<Dim>& mesh,
    const linalg::Vector& solution,
    const problems::Problem& problem,
    const discretization::element::FiniteElement& fe,
    const discretization::quadrature::QuadratureRule& quad,
    const std::vector<const boundary::BoundaryCondition<Dim>*>& bcs,
    const std::string& base_dir = "results")
{
    const fs::path run_dir = fs::path(base_dir) / detail::make_timestamp();
    fs::create_directories(run_dir);

    detail::write_solution_csv<Dim>(run_dir / "solution.csv", mesh, solution);
    detail::write_meta_json<Dim>(run_dir / "meta.json", mesh, problem, fe, quad, bcs);

    return run_dir;
}

} // namespace fem::io