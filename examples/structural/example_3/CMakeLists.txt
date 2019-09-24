add_executable(structural_example_3
                example_3.cpp
                ${MAST_ROOT_DIR}/examples/structural/base/thermal_stress_jacobian_scaling_function.cpp)

target_include_directories(structural_example_3 PRIVATE
                            ${CMAKE_CURRENT_LIST_DIR}
                            ${PROJECT_SOURCE_DIR}
                            ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(structural_example_3 mast)

install(TARGETS structural_example_3
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
