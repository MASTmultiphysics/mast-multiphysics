target_sources(mast_catch_tests
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/mast_boundary_condition_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_thermoelastic_load.cpp)


# BoundaryConditionBase tests that depend on FunctionSetBase tests being successful
add_test(NAME BoundaryConditionBase
    COMMAND mast_catch_tests -w NoTests "boundary_condition_base")
set_tests_properties(BoundaryConditionBase
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED FunctionSetBase
        FIXTURES_SETUP BoundaryConditionBase)

add_test(NAME BoundaryConditionBase_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "boundary_condition_base")
set_tests_properties(BoundaryConditionBase_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED FunctionSetBase_mpi
        FIXTURES_SETUP BoundaryConditionBase_mpi)

# Temperature Boundary Condition Test
add_test(NAME Thermoelastic_Load
    COMMAND mast_catch_tests -w NoTests "create_thermoelastic_load")
set_tests_properties(Thermoelastic_Load
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED BoundaryConditionBase
        FIXTURES_SETUP Thermoelastic_Load)

add_test(NAME Thermoelastic_Load_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "create_thermoelastic_load")
set_tests_properties(Thermoelastic_Load_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED BoundaryConditionBase_mpi
        FIXTURES_SETUP Thermoelastic_Load_mpi)