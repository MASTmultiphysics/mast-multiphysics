# since an optimizer is required for optimization we turn this on only if
# NLOPT was enabled for MAST
if(ENABLE_NLOPT)
    add_executable(structural_example_5
                   example_5.cpp)

    target_include_directories(structural_example_5 PRIVATE
                                ${CMAKE_CURRENT_LIST_DIR}
                                ${PROJECT_SOURCE_DIR}
                                ${PROJECT_SOURCE_DIR}/src)
    target_link_libraries(structural_example_5 mast)

    install(TARGETS structural_example_5
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
endif()
