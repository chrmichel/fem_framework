#pragma once
#include <vector>
#include <cstddef>

namespace fem::linalg {

class SparseMatrix {
public:
    // Aufbau aus vorgegebenem Sparsity-Pattern
    // row_ptr[i]..row_ptr[i+1] = Spaltenindizes für Zeile i
    SparseMatrix(std::size_t n_rows,
                 std::vector<std::size_t> row_ptr,
                 std::vector<std::size_t> col_idx);

    std::size_t rows() const;
    std::size_t cols() const;
    std::size_t nnz()  const;  // Anzahl gespeicherter Einträge

    // Zugriff — wirft wenn (i,j) nicht im Pattern
    double& operator()(std::size_t i, std::size_t j);
    double  operator()(std::size_t i, std::size_t j) const;

    void set_zero();  // alle Werte auf 0, Pattern bleibt

    const std::vector<std::size_t>& row_ptr() const;
    const std::vector<std::size_t>& col_idx() const;
    const std::vector<double>&      values()  const;

    std::vector<double>& values_mutable() { return m_values; }

private:
    std::size_t m_rows;
    std::vector<std::size_t> m_row_ptr;
    std::vector<std::size_t> m_col_idx;
    std::vector<double>      m_values;

    // Hilfsfunktion: findet Position von (i,j) in m_values
    std::size_t find(std::size_t i, std::size_t j) const;
};

} // namespace fem::linalg