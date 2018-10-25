target_sources(mast
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/fe_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/fe_base.h
        ${CMAKE_CURRENT_LIST_DIR}/local_1d_elem.cpp
        ${CMAKE_CURRENT_LIST_DIR}/local_1d_elem.h
        ${CMAKE_CURRENT_LIST_DIR}/local_2d_elem.cpp
        ${CMAKE_CURRENT_LIST_DIR}/local_2d_elem.h
        ${CMAKE_CURRENT_LIST_DIR}/local_3d_elem.cpp
        ${CMAKE_CURRENT_LIST_DIR}/local_3d_elem.h
        ${CMAKE_CURRENT_LIST_DIR}/local_elem_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/local_elem_base.h)

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/mesh
        FILES_MATCHING PATTERN "*.h")
