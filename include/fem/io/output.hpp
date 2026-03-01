#pragma once

#include "fem/core/mesh.hpp"
#include "fem/linalg/vector.hpp"
#include "fem/problems/problem.hpp"
#include "meta.hpp"

#include <filesystem>
#include <functional>
#include <map>
#include <string>
#include <variant>

namespace fem::io {

using ScalarFunction = std::function<double(double)>;

// Creates results/<timestamp>/ and returns the created directory.
std::filesystem::path make_results_dir(const std::filesystem::path& root = "results");

// Writes solution.csv into out_dir. Creates out_dir if needed.
void write_solution_csv(const std::filesystem::path& out_dir,
                        const fem::Mesh& mesh,
                        const fem::linalg::Vector& u,
                        ScalarFunction u_exact = nullptr);

// Create meta from Mesh + Problem + framework defaults
Meta make_meta(const fem::problems::Problem& problem);

// Convenience: write meta.json by only giving problem
void write_meta_json(const std::filesystem::path& out_dir,
                     const fem::problems::Problem& problem);

// Returns newest subdirectory inside root (lexicographically max for our timestamp format).
// Throws if none exists.
std::filesystem::path latest_results_dir(const std::filesystem::path& root = "results");

} // namespace fem::io