add_executable(mesh_test
    mesh_test.cpp
)
target_link_libraries(mesh_test PRIVATE fem)
target_include_directories(mesh_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME mesh_test COMMAND mesh_test)

add_executable(linear_dirichlet_test
    linear_dirichlet_test.cpp
)
target_link_libraries(linear_dirichlet_test PRIVATE fem)
target_include_directories(linear_dirichlet_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME linear_dirichlet_test COMMAND linear_dirichlet_test)

add_executable(l2_error_test
    l2_error_test.cpp
)
target_link_libraries(l2_error_test PRIVATE fem)
target_include_directories(l2_error_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME l2_error_test COMMAND l2_error_test)

add_executable(convergence_test
    convergence_test.cpp
)
target_link_libraries(convergence_test PRIVATE fem)
target_include_directories(convergence_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME convergence_test COMMAND convergence_test)

add_executable(reaction_test
    reaction_test.cpp
)
target_link_libraries(reaction_test PRIVATE fem)
target_include_directories(reaction_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME reaction_test COMMAND reaction_test)

add_executable(output_test
    output_test.cpp
)
target_link_libraries(output_test PRIVATE fem)
target_include_directories(output_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME output_test COMMAND output_test)

add_executable(p2_element_test
    p2_element_test.cpp
)
target_link_libraries(p2_element_test PRIVATE fem)
target_include_directories(p2_element_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME p2_element_test COMMAND p2_element_test)

add_executable(p2_convergence_test
    p2_convergence_test.cpp
)
target_link_libraries(p2_convergence_test PRIVATE fem)
target_include_directories(p2_convergence_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME p2_convergence_test COMMAND p2_convergence_test)
add_executable(neumann_bc_test
    neumann_bc_test.cpp
)
target_link_libraries(neumann_bc_test PRIVATE fem)
target_include_directories(neumann_bc_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME neumann_bc_test COMMAND neumann_bc_test)
add_executable(adaptive_refinement_test
    adaptive_refinement_test.cpp
)
target_link_libraries(adaptive_refinement_test PRIVATE fem)
target_include_directories(adaptive_refinement_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME adaptive_refinement_test COMMAND adaptive_refinement_test)
add_executable(sparsity_test
    sparsity_test.cpp
)
target_link_libraries(sparsity_test PRIVATE fem)
target_include_directories(sparsity_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME sparsity_test COMMAND sparsity_test)
add_executable(sparse_matrix_test
    sparse_matrix_test.cpp
)
target_link_libraries(sparse_matrix_test PRIVATE fem)
target_include_directories(sparse_matrix_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME sparse_matrix_test COMMAND sparse_matrix_test)
add_executable(mesh2d_test
    mesh2d_test.cpp
)
target_link_libraries(mesh2d_test PRIVATE fem)
target_include_directories(mesh2d_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME mesh2d_test COMMAND mesh2d_test)
add_executable(p1_2d_test
    p1_2d_test.cpp
)
target_link_libraries(p1_2d_test PRIVATE fem)
target_include_directories(p1_2d_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME p1_2d_test COMMAND p1_2d_test)
add_executable(triangle_quadrature_test
    triangle_quadrature_test.cpp
)
target_link_libraries(triangle_quadrature_test PRIVATE fem)
target_include_directories(triangle_quadrature_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME triangle_quadrature_test COMMAND triangle_quadrature_test)
add_executable(poisson2d_convergence_test
    poisson2d_convergence_test.cpp
)
target_link_libraries(poisson2d_convergence_test PRIVATE fem)
target_include_directories(poisson2d_convergence_test PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME poisson2d_convergence_test COMMAND poisson2d_convergence_test)