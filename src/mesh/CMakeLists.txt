target_sources(mast
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/fe_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/fe_base.h
        ${CMAKE_CURRENT_LIST_DIR}/geom_elem.cpp
        ${CMAKE_CURRENT_LIST_DIR}/geom_elem.h)

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/mesh
        FILES_MATCHING PATTERN "*.h")
