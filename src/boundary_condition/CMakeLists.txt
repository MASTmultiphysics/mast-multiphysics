target_sources(mast
        PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/dirichlet_boundary_condition.cpp
        ${CMAKE_CURRENT_LIST_DIR}/dirichlet_boundary_condition.h
        ${CMAKE_CURRENT_LIST_DIR}/point_load_condition.cpp
        ${CMAKE_CURRENT_LIST_DIR}/point_load_condition.h)