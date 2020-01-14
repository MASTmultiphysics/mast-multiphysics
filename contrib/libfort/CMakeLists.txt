# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/external
        FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ./ DESTINATION include/external
        FILES_MATCHING PATTERN "*.hpp")