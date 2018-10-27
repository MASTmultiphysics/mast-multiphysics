target_sources(mast
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/basis_matrix_coordinate.cpp
        ${CMAKE_CURRENT_LIST_DIR}/basis_matrix_coordinate.h
        ${CMAKE_CURRENT_LIST_DIR}/coordinate_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/coordinate_base.h
        ${CMAKE_CURRENT_LIST_DIR}/polar_coordinate.cpp
        ${CMAKE_CURRENT_LIST_DIR}/polar_coordinate.h)

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/coordinates
        FILES_MATCHING PATTERN "*.h")