#include "fem/io/output.hpp"

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

namespace fem::io {
namespace fs = std::filesystem;

namespace {

std::string make_timestamp()
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

std::string escape_json(const std::string& s)
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

void write_meta_value(std::ostream& out, const MetaValue& value)
{
    std::visit(
        [&](const auto& v)
        {
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

void write_meta_object(std::ostream& out,
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
        if (count < meta.size()) {
            out << ",";
        }
        out << "\n";
    }

    out << std::string(indent - 2, ' ') << "}";
}

void write_solution_csv(const fs::path& file,
                        const fem::core::Mesh& mesh,
                        const fem::linalg::Vector& solution)
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

    out << "node_id,x,u\n";
    for (std::size_t i = 0; i < mesh.n_nodes(); ++i) {
        out << i << ","
            << std::setprecision(16) << mesh.node(i) << ","
            << std::setprecision(16) << solution(i) << "\n";
    }
}

void write_meta_json(const fs::path& file,
                     const fem::core::Mesh& mesh,
                     const fem::problems::Problem& problem,
                     const fem::discretization::element::FiniteElement& fe,
                     const fem::discretization::quadrature::QuadratureRule& quad,
                     const fem::boundary::BoundaryCondition& bc)
{
    std::ofstream out(file);
    if (!out) {
        throw std::runtime_error(
            "write_meta_json: could not open file " + file.string());
    }

    out << "{\n";
    out << "  \"problem\": \"" << escape_json(problem.name()) << "\",\n";

    out << "  \"problem_meta\": ";
    write_meta_object(out, problem.meta(), 4);
    out << ",\n";

    out << "  \"finite_element\": \"" << escape_json(fe.name()) << "\",\n";
    out << "  \"quadrature\": \"" << escape_json(quad.name()) << "\",\n";
    out << "  \"boundary_condition\": \"" << escape_json(bc.name()) << "\",\n";
    out << "  \"n_nodes\": " << mesh.n_nodes() << ",\n";
    out << "  \"n_elements\": " << mesh.n_elements() << "\n";
    out << "}\n";
}

} // namespace

fs::path write_results(const core::Mesh& mesh,
                       const linalg::Vector& solution,
                       const problems::Problem& problem,
                       const discretization::element::FiniteElement& fe,
                       const discretization::quadrature::QuadratureRule& quad,
                       const boundary::BoundaryCondition& bc,
                       const std::string& base_dir)
{
    const fs::path run_dir = fs::path(base_dir) / make_timestamp();
    fs::create_directories(run_dir);

    write_solution_csv(run_dir / "solution.csv", mesh, solution);
    write_meta_json(run_dir / "meta.json", mesh, problem, fe, quad, bc);

    return run_dir;
}

} // namespace fem::io