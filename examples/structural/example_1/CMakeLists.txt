add_executable(structural_example_1
               example_1.cpp)

target_link_libraries(structural_example_1 mast)

install(TARGETS structural_example_1
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
