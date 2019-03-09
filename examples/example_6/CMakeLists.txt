add_executable(example_6
                example_6.cpp
                ${MAST_ROOT_DIR}/examples/fluid/meshing/naca0012.cpp
                ${MAST_ROOT_DIR}/examples/fluid/meshing/naca0012_wing.cpp
                ${MAST_ROOT_DIR}/examples/fluid/meshing/panel_mesh_2D.cpp
                ${MAST_ROOT_DIR}/examples/fluid/meshing/panel_mesh_3D.cpp
                ${MAST_ROOT_DIR}/examples/fluid/meshing/cylinder.cpp
                ${MAST_ROOT_DIR}/examples/fluid/meshing/mesh_initializer.cpp)

target_include_directories(example_6 PRIVATE
                            ${CMAKE_CURRENT_LIST_DIR}
                            ${PROJECT_SOURCE_DIR}
                            ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(example_6 mast)

install(TARGETS example_6
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)

