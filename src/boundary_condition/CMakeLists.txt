target_sources(mast
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/dirichlet_boundary_condition.cpp
        ${CMAKE_CURRENT_LIST_DIR}/dirichlet_boundary_condition.h
        ${CMAKE_CURRENT_LIST_DIR}/dirichlet_dof_boundary_condition.cpp
        ${CMAKE_CURRENT_LIST_DIR}/dirichlet_dof_boundary_condition.h
        ${CMAKE_CURRENT_LIST_DIR}/point_load_condition.cpp
        ${CMAKE_CURRENT_LIST_DIR}/point_load_condition.h)

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/boundary_condition
        FILES_MATCHING PATTERN "*.h")
