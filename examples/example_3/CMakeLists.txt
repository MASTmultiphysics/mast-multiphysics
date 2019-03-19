add_executable(example_3 
               example_3.cpp 
               constrain_beam_dofs.cpp
               ${MAST_ROOT_DIR}/examples/fluid/meshing/panel_mesh_2D.cpp  
               ${MAST_ROOT_DIR}/examples/fluid/meshing/mesh_initializer.cpp
               ${MAST_ROOT_DIR}/examples/base/augment_ghost_elem_send_list.cpp
               ${MAST_ROOT_DIR}/examples/base/plot_results.cpp)

target_include_directories(example_3 PRIVATE 
                           ${CMAKE_CURRENT_LIST_DIR}
                           ${PROJECT_SOURCE_DIR}
                           ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(example_3 mast)

install(TARGETS example_3
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
