add_executable(example_5
               example_5.cpp)

target_include_directories(example_5 PRIVATE
                            ${CMAKE_CURRENT_LIST_DIR}
                            ${PROJECT_SOURCE_DIR}
                            ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(example_5 mast)

install(TARGETS example_5
RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)

