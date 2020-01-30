if(ENABLE_NASTRANIO)

    add_executable(structural_example_7
            example_7.cpp)

    target_link_libraries(structural_example_7 mast)

    install(TARGETS structural_example_7
            RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)

    configure_file(${CMAKE_CURRENT_LIST_DIR}/example_7_acoss_mesh.bdf ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

    add_test(NAME structural_example_7_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:structural_example_7> -ksp_type "preonly" -pc_type "lu")

endif()
