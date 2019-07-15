add_executable(structural_example_4
                example_4.cpp)

target_link_libraries(structural_example_4 mast)

install(TARGETS structural_example_4
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
