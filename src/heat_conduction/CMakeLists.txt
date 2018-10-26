target_sources(mast
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/heat_conduction_elem_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/heat_conduction_elem_base.h
        ${CMAKE_CURRENT_LIST_DIR}/heat_conduction_nonlinear_assembly.cpp
        ${CMAKE_CURRENT_LIST_DIR}/heat_conduction_nonlinear_assembly.h
        ${CMAKE_CURRENT_LIST_DIR}/heat_conduction_system_initialization.cpp
        ${CMAKE_CURRENT_LIST_DIR}/heat_conduction_system_initialization.h
        ${CMAKE_CURRENT_LIST_DIR}/heat_conduction_transient_assembly.cpp
        ${CMAKE_CURRENT_LIST_DIR}/heat_conduction_transient_assembly.h)

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/heat_conduction
        FILES_MATCHING PATTERN "*.h")