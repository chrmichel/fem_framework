#include "fem/io/output.hpp"
#include "fem/problems/problem.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <ctime>
#include <cstdint>

#ifndef FEM_SOURCE_DIR
#define FEM_SOURCE_DIR "."
#endif 

namespace fem::io {

static void require_size_match(const fem::Mesh& mesh, const fem::linalg::Vector& u)
{
    if (u.size() != mesh.n_nodes()) {
        throw std::invalid_argument("write_solution_csv: u.size() must match mesh.n_nodes().");
    }
}

static std::string timestamp_local()
{
    std::time_t t = std::time(nullptr);
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

std::filesystem::path make_results_dir(const std::filesystem::path& root)
{
    std::filesystem::path base = std::filesystem::path(FEM_SOURCE_DIR) / root;
    std::error_code ec;
    std::filesystem::create_directories(base, ec);
    if (ec) {
        throw std::runtime_error("make_results_dir: cannot create root directory.");
    }

    auto run_dir = base / timestamp_local();
    std::filesystem::create_directories(run_dir, ec);
    if (ec) {
        throw std::runtime_error("make_results_dir: cannot create run directory.");
    }
    return run_dir;
}

void write_solution_csv(const std::filesystem::path& out_dir,
                        const fem::Mesh& mesh,
                        const fem::linalg::Vector& u,
                        ScalarFunction u_exact)
{
    require_size_match(mesh, u);

    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) {
        throw std::runtime_error("write_solution_csv: cannot create output directory.");
    }

    const auto file = out_dir / "solution.csv";
    std::ofstream out(file);
    if (!out) {
        throw std::runtime_error("write_solution_csv: cannot open solution.csv for writing.");
    }

    out << std::setprecision(17);

    if (u_exact) {
        out << "i,x,u,u_exact,error\n";
        for (std::size_t i = 0; i < mesh.n_nodes(); ++i) {
            const double x  = mesh.node(i);
            const double ue = u_exact(x);
            const double err = u(i) - ue;
            out << i << "," << x << "," << u(i) << "," << ue << "," << err << "\n";
        }
    } else {
        out << "i,x,u\n";
        for (std::size_t i = 0; i < mesh.n_nodes(); ++i) {
            out << i << "," << mesh.node(i) << "," << u(i) << "\n";
        }
    }
}

// ---------- JSON helpers ----------

static std::string json_escape(const std::string& s)
{
    std::string out;
    out.reserve(s.size() + 8);

    for (char c : s) {
        switch (c) {
            case '\"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b";  break;
            case '\f': out += "\\f";  break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:
                // keep ASCII control chars out of JSON
                if (static_cast<unsigned char>(c) < 0x20) {
                    std::ostringstream oss;
                    oss << "\\u"
                        << std::hex << std::uppercase << std::setw(4) << std::setfill('0')
                        << (int)(unsigned char)c;
                    out += oss.str();
                } else {
                    out += c;
                }
        }
    }
    return out;
}

static void write_json_value(std::ostream& os, const MetaValue& v)
{
    if (std::holds_alternative<std::string>(v)) {
        os << "\"" << json_escape(std::get<std::string>(v)) << "\"";
    } else if (std::holds_alternative<double>(v)) {
        // ensure consistent formatting
        os << std::setprecision(17) << std::get<double>(v);
    } else if (std::holds_alternative<std::int64_t>(v)) {
        os << std::get<std::int64_t>(v);
    } else if (std::holds_alternative<bool>(v)) {
        os << (std::get<bool>(v) ? "true" : "false");
    }
}

void write_meta_json(const std::filesystem::path& out_dir,
                     const Meta& meta)
{
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) {
        throw std::runtime_error("write_meta_json: cannot create output directory.");
    }

    const auto file = out_dir / "meta.json";
    std::ofstream out(file);
    if (!out) {
        throw std::runtime_error("write_meta_json: cannot open meta.json for writing.");
    }

    out << "{\n";

    bool first = true;
    for (const auto& [k, v] : meta) {
        if (!first) out << ",\n";
        first = false;

        out << "  \"" << json_escape(k) << "\": ";
        write_json_value(out, v);
    }

    out << "\n}\n";
}

Meta make_meta(const fem::problems::Problem& problem)
{
    Meta m;

    const auto& mesh = problem.mesh();

    // framework/general
    m["framework"] = std::string("FEMFramework");
    m["phase"] = std::int64_t(1);
    m["dimension"] = static_cast<std::int64_t>(mesh.dimension());
    m["domain_a"] = mesh.a();
    m["domain_b"] = mesh.b();
    m["elements"] = static_cast<std::int64_t>(mesh.n_elements());
    m["nodes"] = static_cast<std::int64_t>(mesh.n_nodes());

    // discretization/assembly info (Phase 1 fixed)
    m["element_type"] = std::string("P1");
    m["quadrature"] = std::string("Gauss2");
    m["matrix_type"] = std::string("Dense");
    m["solver"] = std::string("GaussianEliminationPivoting");

    // problem-specific
    m["problem"] = problem.name();

    // merge in problem.meta()
    Meta pm = problem.meta();
    for (auto& [k, v] : pm) {
        m[k] = v; // overwrite allowed
    }

    return m;
}

void write_meta_json(const std::filesystem::path& out_dir,
                     const fem::problems::Problem& problem)
{
    write_meta_json(out_dir, make_meta(problem));
}

std::filesystem::path latest_results_dir(const std::filesystem::path& root)
{
    if (!std::filesystem::exists(root) || !std::filesystem::is_directory(root)) {
        throw std::runtime_error("latest_results_dir: results root does not exist or is not a directory.");
    }

    std::filesystem::path best;
    bool found = false;

    for (const auto& entry : std::filesystem::directory_iterator(root)) {
        if (!entry.is_directory()) continue;
        const auto name = entry.path().filename().string();

        // Our timestamp format YYYY-MM-DD_HH-MM-SS is lexicographically sortable
        if (!found || name > best.filename().string()) {
            best = entry.path();
            found = true;
        }
    }

    if (!found) {
        throw std::runtime_error("latest_results_dir: no result directories found.");
    }
    return best;
}


} // namespace fem::io