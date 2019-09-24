add_executable(structural_example_2
               example_2.cpp)

target_include_directories(structural_example_2 PRIVATE
                            ${CMAKE_CURRENT_LIST_DIR}
                            ${PROJECT_SOURCE_DIR}
                            ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(structural_example_2 mast)

install(TARGETS structural_example_2
RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
