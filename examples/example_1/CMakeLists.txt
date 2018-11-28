add_executable(example_1
               ex1_axial_bar_extension.cpp)

target_link_libraries(example_1 mast)

install(TARGETS example_1
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
