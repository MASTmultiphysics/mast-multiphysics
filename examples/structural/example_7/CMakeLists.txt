if(ENABLE_NASTRANIO)

    add_executable(structural_example_7
            example_7.cpp)

    target_link_libraries(structural_example_7 mast)

# Move BDF to binary directory for convenient testing.
configure_file(${CMAKE_CURRENT_LIST_DIR}/example_7_acoss_mesh.bdf ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

## Test on single processor, default options.
#add_test(NAME structural_example_7
#        COMMAND $<TARGET_FILE:structural_example_7>)
#
## Test multiple processors, parallel direct solver using external MUMPS package.
#add_test(NAME structural_example_7_mpi
#    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:structural_example_7>
#            -st_ksp_type gmres -st_pc_type bjacobi -st_ksp_rtol 1e-9)

endif()