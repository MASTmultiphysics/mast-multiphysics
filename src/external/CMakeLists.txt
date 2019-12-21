target_sources(mast
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/fort.c
        ${CMAKE_CURRENT_LIST_DIR}/fort.h
        ${CMAKE_CURRENT_LIST_DIR}/fort.hpp)

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/external
        FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ./ DESTINATION include/external
        FILES_MATCHING PATTERN "*.hpp")