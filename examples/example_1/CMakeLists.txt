add_executable(ex1_axial_bar_extension
    ex1_axial_bar_extension.cpp)

target_link_libraries(ex1_axial_bar_extension mast)

install(TARGETS ex1_axial_bar_extension
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)