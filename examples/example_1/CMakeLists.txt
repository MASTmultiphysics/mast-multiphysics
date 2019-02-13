add_executable(example_1
               example_1.cpp)

target_link_libraries(example_1 mast)

install(TARGETS example_1
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
