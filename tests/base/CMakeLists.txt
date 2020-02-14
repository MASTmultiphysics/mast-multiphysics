target_sources(mast_catch_tests
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/mast_parameter.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_constant_field_function.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_function_set_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_mesh.cpp)

# FIXME: MPI tests seem to either run very slow or hang up intermittently
# This has occured in:
#   libMesh_Mesh_Generation_1d_mpi  - hung up
        
# libMesh mesh generation test
add_test(NAME libMesh_Mesh_Generation_1d
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "libmesh_mesh_generation_1d")
set_tests_properties(libMesh_Mesh_Generation_1d
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_SETUP libMesh_Mesh_Generation_1d)
        
add_test(NAME libMesh_Mesh_Generation_1d_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "libmesh_mesh_generation_1d")
set_tests_properties(libMesh_Mesh_Generation_1d_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_SETUP libMesh_Mesh_Generation_1d_mpi)

add_test(NAME libMesh_Mesh_Generation_2d
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "libmesh_mesh_generation_2d")
set_tests_properties(libMesh_Mesh_Generation_2d
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_SETUP libMesh_Mesh_Generation_2d)
        
add_test(NAME libMesh_Mesh_Generation_2d_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "libmesh_mesh_generation_2d")
set_tests_properties(libMesh_Mesh_Generation_2d_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_SETUP libMesh_Mesh_Generation_2d_mpi)

# Parameter tests
add_test(NAME Parameter
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests parameters)
set_tests_properties(Parameter
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_SETUP Parameter)
    
add_test(NAME Parameter_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests parameters)
set_tests_properties(Parameter_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_SETUP Parameter_mpi)

# Constant Field Function tests that depend on parameter tests being successful
add_test(NAME ConstantFieldFunction 
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_field_functions)
set_tests_properties(ConstantFieldFunction 
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Parameter
        FIXTURES_SETUP ConstantFieldFunction)
        
add_test(NAME ConstantFieldFunction_mpi 
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_field_functions)
set_tests_properties(ConstantFieldFunction_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Parameter_mpi
        FIXTURES_SETUP ConstantFieldFunction_mpi)
                     
# FunctionSetBase tests that depend on constant field function tests being successful
add_test(NAME FunctionSetBase
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "function_set_base")
set_tests_properties(FunctionSetBase
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED ConstantFieldFunction
        FIXTURES_SETUP FunctionSetBase)

add_test(NAME FunctionSetBase_mpi
COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "function_set_base")
set_tests_properties(FunctionSetBase_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED ConstantFieldFunction_mpi
        FIXTURES_SETUP FunctionSetBase_mpi)
