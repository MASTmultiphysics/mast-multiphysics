add_executable(structural_example_9
    example_9.cpp)

target_link_libraries(structural_example_9 mast)

install(TARGETS structural_example_9
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)

# Test on single processor, PETSc built-in LU direct linear solver (sequential).
add_test(NAME structural_example_9
    COMMAND $<TARGET_FILE:structural_example_9> -ksp_type preonly -pc_type lu -options_view)
set_tests_properties(structural_example_9
    PROPERTIES
        LABELS "SEQ")

# Test multiple processors, parallel direct solver using external MUMPS package.
add_test(NAME structural_example_9_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:structural_example_9>
            -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -options_view)
set_tests_properties(structural_example_9_mpi
    PROPERTIES
        LABELS "MPI")
