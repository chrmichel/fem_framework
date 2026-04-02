add_library(fem)

target_sources(fem PRIVATE
    core/mesh.cpp
    core/mesh_utils.cpp
    assembly/assembler.cpp
    boundary/dirichlet_bc.cpp
    boundary/neumann_bc.cpp
    linalg/matrix.cpp
    linalg/vector.cpp
    linalg/solver.cpp
    solve/driver.cpp
    io/output.cpp
    discretization/element/lagrange_p1_1d.cpp
    discretization/element/lagrange_p2_1d.cpp
    discretization/quadrature/gauss_legendre_1d.cpp
    problems/poisson_1d.cpp
)

target_include_directories(fem
    PUBLIC
        ${CMAKE_SOURCE_DIR}/include
)

target_compile_features(fem PUBLIC cxx_std_20)
target_compile_definitions(fem PUBLIC
    FEM_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)