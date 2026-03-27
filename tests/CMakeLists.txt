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