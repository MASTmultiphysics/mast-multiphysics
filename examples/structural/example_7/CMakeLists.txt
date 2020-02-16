if(ENABLE_NASTRANIO)

    add_executable(structural_example_7
        example_7.cpp)

    target_link_libraries(structural_example_7 mast)

    # Move BDF to binary directory for convenient testing.
    configure_file(${CMAKE_CURRENT_LIST_DIR}/example_7_acoss_mesh.bdf ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

    install(TARGETS structural_example_7
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)

    # Test on single processor, PETSc built-in LU direct linear solver (sequential).
    add_test(NAME structural_example_7
        COMMAND $<TARGET_FILE:structural_example_7> -ksp_type preonly -pc_type lu -options_view)
    set_tests_properties(structural_example_7
        PROPERTIES
            LABELS "SHORT;SEQ")

    # Test multiple processors, parallel direct solver using external MUMPS package.
    add_test(NAME structural_example_7_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:structural_example_7>
                -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -options_view)
    set_tests_properties(structural_example_7_mpi
        PROPERTIES
            LABELS "SHORT;MPI")

endif()