target_sources( mast_catch_tests
                PRIVATE
                ${CMAKE_CURRENT_LIST_DIR}/mast_parameter.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_constant_field_function.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_function_set_base.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_mesh.cpp)

# libMesh mesh generation test
add_test(NAME libMesh_Mesh_Generation_1d
         COMMAND mast_catch_tests "libmesh_mesh_generation_1d")
set_tests_properties(libMesh_Mesh_Generation_1d
                     PROPERTIES
                     FIXTURES_SETUP libMesh_Mesh_Generation_1d)

add_test(NAME libMesh_Mesh_Generation_2d
         COMMAND mast_catch_tests "libmesh_mesh_generation_2d")
set_tests_properties(libMesh_Mesh_Generation_2d
                     PROPERTIES
                     FIXTURES_SETUP libMesh_Mesh_Generation_2d)

# Parameter tests
add_test(NAME Parameter 
         COMMAND mast_catch_tests parameters)
set_tests_properties(Parameter PROPERTIES FIXTURES_SETUP Parameter)

# Constant Field Function tests that depend on parameter tests being successful
add_test(NAME ConstantFieldFunction 
         COMMAND mast_catch_tests constant_field_functions)
set_tests_properties(ConstantFieldFunction 
                     PROPERTIES 
                     FIXTURES_REQUIRED Parameter
                     FIXTURES_SETUP ConstantFieldFunction)
                     
# FunctionSetBase tests that depend on constant field function tests being successful
add_test(NAME FunctionSetBase
         COMMAND mast_catch_tests "function_set_base")
set_tests_properties(FunctionSetBase
                     PROPERTIES 
                     FIXTURES_REQUIRED ConstantFieldFunction
                     FIXTURES_SETUP FunctionSetBase)


