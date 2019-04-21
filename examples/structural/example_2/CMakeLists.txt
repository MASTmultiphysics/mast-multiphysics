add_executable(structural_example_2
               example_2.cpp)

target_link_libraries(structural_example_2 mast)

install(TARGETS structural_example_2
RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
