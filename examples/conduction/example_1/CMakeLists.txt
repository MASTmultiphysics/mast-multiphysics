add_executable(conduction_example_1
                example_1.cpp)

target_include_directories(conduction_example_1 PRIVATE
                            ${CMAKE_CURRENT_LIST_DIR}
                            ${PROJECT_SOURCE_DIR}
                            ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(conduction_example_1 mast)

install(TARGETS conduction_example_1
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)

