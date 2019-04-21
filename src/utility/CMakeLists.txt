target_sources(mast
                PRIVATE
                ${CMAKE_CURRENT_LIST_DIR}/plot.cpp
                ${CMAKE_CURRENT_LIST_DIR}/plot.h)

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/solver
        FILES_MATCHING PATTERN "*.h")
