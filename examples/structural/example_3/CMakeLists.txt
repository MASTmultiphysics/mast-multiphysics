add_executable(structural_example_3
                example_3.cpp)

target_link_libraries(structural_example_3 mast)

install(TARGETS structural_example_3
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
