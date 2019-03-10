add_executable(fluid_example_1
                example_1.cpp
                ${MAST_ROOT_DIR}/examples/fluid/meshing/naca0012.cpp
                ${MAST_ROOT_DIR}/examples/fluid/meshing/naca0012_wing.cpp
                ${MAST_ROOT_DIR}/examples/fluid/meshing/panel_mesh_2D.cpp
                ${MAST_ROOT_DIR}/examples/fluid/meshing/panel_mesh_3D.cpp
                ${MAST_ROOT_DIR}/examples/fluid/meshing/cylinder.cpp
                ${MAST_ROOT_DIR}/examples/fluid/meshing/mesh_initializer.cpp)

target_include_directories(fluid_example_1 PRIVATE
                            ${CMAKE_CURRENT_LIST_DIR}
                            ${PROJECT_SOURCE_DIR}
                            ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(fluid_example_1 mast)

install(TARGETS fluid_example_1
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)

