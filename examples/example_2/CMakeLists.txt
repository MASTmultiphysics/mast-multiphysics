add_executable(example_2
               ex2_plate_static.cpp)

target_link_libraries(example_2 mast)

install(TARGETS example_2
RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
