target_sources(mast
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/assembly_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/assembly_base.h
        ${CMAKE_CURRENT_LIST_DIR}/assembly_elem_operation.cpp
        ${CMAKE_CURRENT_LIST_DIR}/assembly_elem_operation.h
        ${CMAKE_CURRENT_LIST_DIR}/boundary_condition_base.h
        ${CMAKE_CURRENT_LIST_DIR}/complex_assembly_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/complex_assembly_base.h
        ${CMAKE_CURRENT_LIST_DIR}/complex_assembly_elem_operations.cpp
        ${CMAKE_CURRENT_LIST_DIR}/complex_assembly_elem_operations.h
        ${CMAKE_CURRENT_LIST_DIR}/complex_mesh_field_function.cpp
        ${CMAKE_CURRENT_LIST_DIR}/complex_mesh_field_function.h
        ${CMAKE_CURRENT_LIST_DIR}/constant_field_function.cpp
        ${CMAKE_CURRENT_LIST_DIR}/constant_field_function.h
        ${CMAKE_CURRENT_LIST_DIR}/dof_coupling_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/dof_coupling_base.h
        ${CMAKE_CURRENT_LIST_DIR}/eigenproblem_assembly.cpp
        ${CMAKE_CURRENT_LIST_DIR}/eigenproblem_assembly.h
        ${CMAKE_CURRENT_LIST_DIR}/eigenproblem_assembly_elem_operations.cpp
        ${CMAKE_CURRENT_LIST_DIR}/eigenproblem_assembly_elem_operations.h
        ${CMAKE_CURRENT_LIST_DIR}/elem_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/elem_base.h
        ${CMAKE_CURRENT_LIST_DIR}/field_function_base.h
        ${CMAKE_CURRENT_LIST_DIR}/function_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/function_base.h
        ${CMAKE_CURRENT_LIST_DIR}/function_set_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/function_set_base.h
        ${CMAKE_CURRENT_LIST_DIR}/mast_data_types.h
        ${CMAKE_CURRENT_LIST_DIR}/mesh_field_function.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mesh_field_function.h
        ${CMAKE_CURRENT_LIST_DIR}/nonlinear_implicit_assembly.cpp
        ${CMAKE_CURRENT_LIST_DIR}/nonlinear_implicit_assembly.h
        ${CMAKE_CURRENT_LIST_DIR}/nonlinear_implicit_assembly_elem_operations.cpp
        ${CMAKE_CURRENT_LIST_DIR}/nonlinear_implicit_assembly_elem_operations.h
        ${CMAKE_CURRENT_LIST_DIR}/nonlinear_system.cpp
        ${CMAKE_CURRENT_LIST_DIR}/nonlinear_system.h
        ${CMAKE_CURRENT_LIST_DIR}/output_assembly_elem_operations.cpp
        ${CMAKE_CURRENT_LIST_DIR}/output_assembly_elem_operations.h
        ${CMAKE_CURRENT_LIST_DIR}/parameter.h
        ${CMAKE_CURRENT_LIST_DIR}/physics_discipline_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/physics_discipline_base.h
        ${CMAKE_CURRENT_LIST_DIR}/system_initialization.cpp
        ${CMAKE_CURRENT_LIST_DIR}/system_initialization.h
        ${CMAKE_CURRENT_LIST_DIR}/transient_assembly.cpp
        ${CMAKE_CURRENT_LIST_DIR}/transient_assembly.h
        ${CMAKE_CURRENT_LIST_DIR}/transient_assembly_elem_operations.cpp
        ${CMAKE_CURRENT_LIST_DIR}/transient_assembly_elem_operations.h
        ${CMAKE_CURRENT_LIST_DIR}/warping_assembly.cpp
        ${CMAKE_CURRENT_LIST_DIR}/warping_assembly.h)

configure_file(${CMAKE_CURRENT_LIST_DIR}/mast_config.h.in
               ${CMAKE_CURRENT_LIST_DIR}/mast_config.h)

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/base
        FILES_MATCHING PATTERN "*.h")
