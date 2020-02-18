target_sources( mast_catch_tests
                PRIVATE
                ${CMAKE_CURRENT_LIST_DIR}/mast_thermoelastic_load.cpp)


# # BoundaryConditionBase tests that depend on FunctionSetBase tests being successful
# add_test(NAME BoundaryConditionBase
#          COMMAND mast_catch_tests "boundary_condition_base")
# set_tests_properties(BoundaryConditionBase
#                      PROPERTIES 
#                      FIXTURES_REQUIRED FunctionSetBase
#                      FIXTURES_SETUP BoundaryConditionBase)


# Temperature Boundary Condition Test
add_test(NAME Thermoelastic_Load
         COMMAND mast_catch_tests "create_thermoelastic_load")
set_tests_properties(Thermoelastic_Load
                     PROPERTIES 
                     FIXTURES_REQUIRED BoundaryConditionBase
                     FIXTURES_SETUP Thermoelastic_Load)
