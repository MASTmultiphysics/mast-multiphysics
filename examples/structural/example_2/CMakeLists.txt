add_executable(structural_example_2
    example_2.cpp)

target_link_libraries(structural_example_2 mast)

target_include_directories(structural_example_2
        PRIVATE
        ${CMAKE_SOURCE_DIR})

install(TARGETS structural_example_2
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)

# Test on single processor, PETSc built-in LU direct linear solver (sequential).
add_test(NAME structural_example_2
    COMMAND $<TARGET_FILE:structural_example_2> -ksp_type preonly -pc_type lu -options_view)
set_tests_properties(structural_example_2
    PROPERTIES
        LABELS "SHORT;SEQ")

# Test multiple processors, parallel direct solver using external MUMPS package.
add_test(NAME structural_example_2_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:structural_example_2>
            -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -options_view)
set_tests_properties(structural_example_2_mpi
    PROPERTIES
        LABELS "SHORT;MPI")
