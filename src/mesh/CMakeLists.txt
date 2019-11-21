target_sources(mast
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/fe_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/fe_base.h
        ${CMAKE_CURRENT_LIST_DIR}/geom_elem.cpp
        ${CMAKE_CURRENT_LIST_DIR}/geom_elem.h
        ${CMAKE_CURRENT_LIST_DIR}/mesh_coupling_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mesh_coupling_base.h
        ${CMAKE_CURRENT_LIST_DIR}/pynastranIO.h
        ${CMAKE_CURRENT_LIST_DIR}/pynastranIO.cpp
        ${CMAKE_CURRENT_LIST_DIR}/NastranIO.cpp
        ${CMAKE_CURRENT_LIST_DIR}/NastranIO.h)

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/mesh
        FILES_MATCHING PATTERN "*.h")
