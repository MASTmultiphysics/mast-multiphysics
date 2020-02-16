add_executable(structural_example_3
    example_3.cpp
    ${CMAKE_SOURCE_DIR}/examples/structural/base/thermal_stress_jacobian_scaling_function.cpp)

target_link_libraries(structural_example_3 mast)

target_include_directories(structural_example_3
    PRIVATE
        ${CMAKE_SOURCE_DIR})

install(TARGETS structural_example_3
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)

# Test on single processor, PETSc built-in LU direct linear solver (sequential).
add_test(NAME structural_example_3
        COMMAND $<TARGET_FILE:structural_example_3> -ksp_type preonly -pc_type lu -options_view)
set_tests_properties(structural_example_3
        PROPERTIES
        LABELS "SEQ")

# Test multiple processors, parallel direct solver using external MUMPS package.
add_test(NAME structural_example_3_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:structural_example_3>
        -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -options_view)
set_tests_properties(structural_example_3_mpi
        PROPERTIES
        LABELS "MPI")
