target_sources( mast_catch_tests
                PRIVATE
                ${CMAKE_CURRENT_LIST_DIR}/mast_isotropic_material.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_orthotropic_material.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_isotropic_material_sensitivity.cpp)

# FIXME: MPI tests seem to either run very slow or hang up intermittently
# This has occured in:
#   FunctionSetBase_mpi - Slow
#   Isotropic_Material_2D_Thermoelastic_mpi - Hung Up 
#   Isotropic_Material_3D_Thermoelastic_mpi - Hung Up
                
# ============================================================================
#                          ISOTROPIC MATERIAL TESTS
# ============================================================================
                
# 1D Isotropic Material Tests
#-----------------------------------------------------------------------------
add_test(NAME Isotropic_Material_1D_Thermoelastic
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_thermoelastic_material_1d)
set_tests_properties(Isotropic_Material_1D_Thermoelastic
                    PROPERTIES 
                    LABELS "SEQ"
                    FIXTURES_REQUIRED   ConstantFieldFunction
                    FIXTURES_SETUP      Isotropic_Material_1D_Thermoelastic)
                    
add_test(NAME Isotropic_Material_1D_Thermoelastic_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_thermoelastic_material_1d)
set_tests_properties(Isotropic_Material_1D_Thermoelastic_mpi
                    PROPERTIES 
                    LABELS "MPI"
                    FIXTURES_REQUIRED   ConstantFieldFunction_mpi
                    FIXTURES_SETUP      Isotropic_Material_1D_Thermoelastic_mpi)

add_test(NAME Isotropic_Material_1D_Structural
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_structural_material_1d)
set_tests_properties(Isotropic_Material_1D_Structural
                    PROPERTIES 
                    LABELS "SEQ"
                    FIXTURES_REQUIRED   ConstantFieldFunction
                    FIXTURES_SETUP      Isotropic_Material_1D_Structural)
                    
add_test(NAME Isotropic_Material_1D_Structural_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_structural_material_1d)
set_tests_properties(Isotropic_Material_1D_Structural_mpi
                    PROPERTIES 
                    LABELS "MPI"
                    FIXTURES_REQUIRED   ConstantFieldFunction_mpi
                    FIXTURES_SETUP      Isotropic_Material_1D_Structural_mpi)

add_test(NAME Isotropic_Material_1D_Heat_Transfer
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_heat_transfer_material_1d)
set_tests_properties(Isotropic_Material_1D_Heat_Transfer
                    PROPERTIES 
                    LABELS "SEQ"
                    FIXTURES_REQUIRED   ConstantFieldFunction
                    FIXTURES_SETUP      Isotropic_Material_1D_Heat_Transfer)
                    
add_test(NAME Isotropic_Material_1D_Heat_Transfer_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_heat_transfer_material_1d)
set_tests_properties(Isotropic_Material_1D_Heat_Transfer_mpi
                    PROPERTIES 
                    LABELS "MPI"
                    FIXTURES_REQUIRED   ConstantFieldFunction_mpi
                    FIXTURES_SETUP      Isotropic_Material_1D_Heat_Transfer_mpi)
                    
add_test(NAME Isotropic_Material_1D_Transient_Heat_Transfer
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_transient_heat_transfer_material_1d)
set_tests_properties(Isotropic_Material_1D_Transient_Heat_Transfer
                    PROPERTIES 
                    LABELS "SEQ"
                    FIXTURES_REQUIRED   ConstantFieldFunction
                    FIXTURES_SETUP      Isotropic_Material_1D_Transient_Heat_Transfer)
                    
add_test(NAME Isotropic_Material_1D_Transient_Heat_Transfer_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_transient_heat_transfer_material_1d)
set_tests_properties(Isotropic_Material_1D_Transient_Heat_Transfer_mpi
                    PROPERTIES 
                    LABELS "MPI"
                    FIXTURES_REQUIRED   ConstantFieldFunction_mpi
                    FIXTURES_SETUP      Isotropic_Material_1D_Transient_Heat_Transfer_mpi)
                    
                    
# 2D Isotropic Material Tests
#-----------------------------------------------------------------------------
add_test(NAME Isotropic_Material_2D_Thermoelastic
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_thermoelastic_material_2d)
set_tests_properties(Isotropic_Material_2D_Thermoelastic
                    PROPERTIES 
                    LABELS "SEQ"
                    FIXTURES_REQUIRED   ConstantFieldFunction
                    FIXTURES_SETUP      Isotropic_Material_2D_Thermoelastic)

add_test(NAME Isotropic_Material_2D_Thermoelastic_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_thermoelastic_material_2d)
set_tests_properties(Isotropic_Material_2D_Thermoelastic_mpi
                    PROPERTIES 
                    LABELS "MPI"
                    FIXTURES_REQUIRED   ConstantFieldFunction_mpi
                    FIXTURES_SETUP      Isotropic_Material_2D_Thermoelastic_mpi)
                    
add_test(NAME Isotropic_Material_2D_Structural
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_structural_material_2d)
set_tests_properties(Isotropic_Material_2D_Structural
                    PROPERTIES 
                    LABELS "SEQ"
                    FIXTURES_REQUIRED   ConstantFieldFunction
                    FIXTURES_SETUP      Isotropic_Material_2D_Structural)
                    
add_test(NAME Isotropic_Material_2D_Structural_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_structural_material_2d)
set_tests_properties(Isotropic_Material_2D_Structural_mpi
                    PROPERTIES 
                    LABELS "MPI"
                    FIXTURES_REQUIRED   ConstantFieldFunction_mpi
                    FIXTURES_SETUP      Isotropic_Material_2D_Structural_mpi)

add_test(NAME Isotropic_Material_2D_Heat_Transfer
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_heat_transfer_material_2d)
set_tests_properties(Isotropic_Material_2D_Heat_Transfer
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED   ConstantFieldFunction
                     FIXTURES_SETUP      Isotropic_Material_2D_Heat_Transfer)
                     
add_test(NAME Isotropic_Material_2D_Heat_Transfer_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_heat_transfer_material_2d)
set_tests_properties(Isotropic_Material_2D_Heat_Transfer_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED   ConstantFieldFunction_mpi
                     FIXTURES_SETUP      Isotropic_Material_2D_Heat_Transfer_mpi)
                    
add_test(NAME Isotropic_Material_2D_Transient_Heat_Transfer
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_transient_heat_transfer_material_2d)
set_tests_properties(Isotropic_Material_2D_Transient_Heat_Transfer
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED   ConstantFieldFunction
                     FIXTURES_SETUP      Isotropic_Material_2D_Transient_Heat_Transfer)
                    
add_test(NAME Isotropic_Material_2D_Transient_Heat_Transfer_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_transient_heat_transfer_material_2d)
set_tests_properties(Isotropic_Material_2D_Transient_Heat_Transfer_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED   ConstantFieldFunction_mpi
                     FIXTURES_SETUP      Isotropic_Material_2D_Transient_Heat_Transfer_mpi)
                    
# 3D Isotropic Material Tests
#-----------------------------------------------------------------------------
add_test(NAME Isotropic_Material_3D_Thermoelastic
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_thermoelastic_material_3d)
set_tests_properties(Isotropic_Material_3D_Thermoelastic
                    PROPERTIES 
                    LABELS "SEQ"
                    FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Isotropic_Material_3D_Thermoelastic_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_thermoelastic_material_3d)
set_tests_properties(Isotropic_Material_3D_Thermoelastic_mpi
                    PROPERTIES 
                    LABELS "MPI"
                    FIXTURES_REQUIRED   ConstantFieldFunction_mpi)
                    
add_test(NAME Isotropic_Material_3D_Structural
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_structural_material_3d)
set_tests_properties(Isotropic_Material_3D_Structural
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED  ConstantFieldFunction)

add_test(NAME Isotropic_Material_3D_Structural_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_structural_material_3d)
set_tests_properties(Isotropic_Material_3D_Structural_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED  ConstantFieldFunction_mpi)
                     
add_test(NAME Isotropic_Material_3D_Heat_Transfer
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_heat_transfer_material_3d)
set_tests_properties(Isotropic_Material_3D_Heat_Transfer
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED  ConstantFieldFunction)
                     
add_test(NAME Isotropic_Material_3D_Heat_Transfer_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_heat_transfer_material_3d)
set_tests_properties(Isotropic_Material_3D_Heat_Transfer_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED  ConstantFieldFunction_mpi)
                     
add_test(NAME Isotropic_Material_3D_Transient_Heat_Transfer
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_transient_heat_transfer_material_3d)
set_tests_properties(Isotropic_Material_3D_Transient_Heat_Transfer
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Isotropic_Material_3D_Transient_Heat_Transfer_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_transient_heat_transfer_material_3d)
set_tests_properties(Isotropic_Material_3D_Transient_Heat_Transfer_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED   ConstantFieldFunction_mpi)
                     
add_test(NAME Isotropic_Material_3D_Dynamic
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_dynamic_material_3d)
set_tests_properties(Isotropic_Material_3D_Dynamic
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED  ConstantFieldFunction)
                     
add_test(NAME Isotropic_Material_3D_Dynamic_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_isotropic_dynamic_material_3d)
set_tests_properties(Isotropic_Material_3D_Dynamic_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED  ConstantFieldFunction_mpi)
                    
                    
# ============================================================================
#                          ORTHOTROPIC MATERIAL TESTS
# ============================================================================

# 1D Orthotropic Material Tests
#-----------------------------------------------------------------------------
add_test(NAME Orthotropic_Material_1D_Thermoelastic
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_thermoelastic_material_1d)
set_tests_properties(Orthotropic_Material_1D_Thermoelastic
                    PROPERTIES 
                    LABELS "SEQ"
                    FIXTURES_REQUIRED   ConstantFieldFunction)
                    
add_test(NAME Orthotropic_Material_1D_Thermoelastic_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_thermoelastic_material_1d)
set_tests_properties(Orthotropic_Material_1D_Thermoelastic_mpi
                    PROPERTIES 
                    LABELS "MPI"
                    FIXTURES_REQUIRED   ConstantFieldFunction_mpi)
                    
add_test(NAME Orthotropic_Material_1D_Structural
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_structural_material_1d)
set_tests_properties(Orthotropic_Material_1D_Structural
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Orthotropic_Material_1D_Structural_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_structural_material_1d)
set_tests_properties(Orthotropic_Material_1D_Structural_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED   ConstantFieldFunction_mpi)
                     
add_test(NAME Orthotropic_Material_1D_Heat_Transfer
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_heat_transfer_material_1d)
set_tests_properties(Orthotropic_Material_1D_Heat_Transfer
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED   ConstantFieldFunction)
                     
add_test(NAME Orthotropic_Material_1D_Heat_Transfer_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_heat_transfer_material_1d)
set_tests_properties(Orthotropic_Material_1D_Heat_Transfer_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED   ConstantFieldFunction_mpi)
                    
add_test(NAME Orthotropic_Material_1D_Transient_Heat_Transfer
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_transient_heat_transfer_material_1d)
set_tests_properties(Orthotropic_Material_1D_Transient_Heat_Transfer
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Orthotropic_Material_1D_Transient_Heat_Transfer_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_transient_heat_transfer_material_1d)
set_tests_properties(Orthotropic_Material_1D_Transient_Heat_Transfer_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED   ConstantFieldFunction_mpi)
                     
# 2D Orthotropic Material Tests
#-----------------------------------------------------------------------------
add_test(NAME Orthotropic_Material_2D_Thermoelastic
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_thermoelastic_material_2d)
set_tests_properties(Orthotropic_Material_2D_Thermoelastic
                    PROPERTIES 
                    LABELS "SEQ"
                    FIXTURES_REQUIRED   ConstantFieldFunction)
                    
add_test(NAME Orthotropic_Material_2D_Thermoelastic_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_thermoelastic_material_2d)
set_tests_properties(Orthotropic_Material_2D_Thermoelastic_mpi
                    PROPERTIES 
                    LABELS "MPI"
                    FIXTURES_REQUIRED   ConstantFieldFunction_mpi)
                    
add_test(NAME Orthotropic_Material_2D_Structural
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_structural_material_2d)
set_tests_properties(Orthotropic_Material_2D_Structural
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Orthotropic_Material_2D_Structural_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_structural_material_2d)
set_tests_properties(Orthotropic_Material_2D_Structural_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED   ConstantFieldFunction_mpi)
                     
add_test(NAME Orthotropic_Material_2D_Heat_Transfer
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_heat_transfer_material_2d)
set_tests_properties(Orthotropic_Material_2D_Heat_Transfer
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Orthotropic_Material_2D_Heat_Transfer_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_heat_transfer_material_2d)
set_tests_properties(Orthotropic_Material_2D_Heat_Transfer_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED   ConstantFieldFunction_mpi)
                    
add_test(NAME Orthotropic_Material_2D_Transient_Heat_Transfer
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_transient_heat_transfer_material_2d)
set_tests_properties(Orthotropic_Material_2D_Transient_Heat_Transfer
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Orthotropic_Material_2D_Transient_Heat_Transfer_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_transient_heat_transfer_material_2d)
set_tests_properties(Orthotropic_Material_2D_Transient_Heat_Transfer_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED   ConstantFieldFunction_mpi)
                     
                     
# 3D Orthotropic Material Tests
#-----------------------------------------------------------------------------
add_test(NAME Orthotropic_Material_3D_Thermoelastic
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_thermoelastic_material_3d)
set_tests_properties(Orthotropic_Material_3D_Thermoelastic
                    PROPERTIES 
                    LABELS "SEQ"
                    FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Orthotropic_Material_3D_Thermoelastic_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_thermoelastic_material_3d)
set_tests_properties(Orthotropic_Material_3D_Thermoelastic_mpi
                    PROPERTIES 
                    LABELS "MPI"
                    FIXTURES_REQUIRED   ConstantFieldFunction_mpi)
                    
add_test(NAME Orthotropic_Material_3D_Structural
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_structural_material_3d)
set_tests_properties(Orthotropic_Material_3D_Structural
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED  ConstantFieldFunction)

add_test(NAME Orthotropic_Material_3D_Structural_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_structural_material_3d)
set_tests_properties(Orthotropic_Material_3D_Structural_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED  ConstantFieldFunction_mpi)
                     
add_test(NAME Orthotropic_Material_3D_Heat_Transfer
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_heat_transfer_material_3d)
set_tests_properties(Orthotropic_Material_3D_Heat_Transfer
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED  ConstantFieldFunction)
                     
add_test(NAME Orthotropic_Material_3D_Heat_Transfer_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_heat_transfer_material_3d)
set_tests_properties(Orthotropic_Material_3D_Heat_Transfer_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED  ConstantFieldFunction_mpi)
                     
add_test(NAME Orthotropic_Material_3D_Transient_Heat_Transfer
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_transient_heat_transfer_material_3d)
set_tests_properties(Orthotropic_Material_3D_Transient_Heat_Transfer
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Orthotropic_Material_3D_Transient_Heat_Transfer_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_transient_heat_transfer_material_3d)
set_tests_properties(Orthotropic_Material_3D_Transient_Heat_Transfer_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED   ConstantFieldFunction_mpi)
                     
add_test(NAME Orthotropic_Material_3D_Dynamic
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_dynamic_material_3d)
set_tests_properties(Orthotropic_Material_3D_Dynamic
                     PROPERTIES 
                     LABELS "SEQ"
                     FIXTURES_REQUIRED  ConstantFieldFunction)

add_test(NAME Orthotropic_Material_3D_Dynamic_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests constant_orthotropic_dynamic_material_3d)
set_tests_properties(Orthotropic_Material_3D_Dynamic_mpi
                     PROPERTIES 
                     LABELS "MPI"
                     FIXTURES_REQUIRED  ConstantFieldFunction_mpi)
    
         
# # ============================================================================
# #                          ANISOTROPIC MATERIAL TESTS
# # ============================================================================
# 
# # 1D Anisotropic Material Tests
# #-----------------------------------------------------------------------------
# add_test(NAME Anisotropic_Material_1D_Thermoelastic
#          COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_anisotropic_thermoelastic_material_1d)
# set_tests_properties(Anisotropic_Material_1D_Thermoelastic
#                     PROPERTIES 
#                     FIXTURES_REQUIRED   ConstantFieldFunction)
#                     
# add_test(NAME Anisotropic_Material_1D_Structural
#          COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_anisotropic_structural_material_1d)
# set_tests_properties(Anisotropic_Material_1D_Structural
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
# 
# add_test(NAME Anisotropic_Material_1D_Heat_Transfer
#          COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_anisotropic_heat_transfer_material_1d)
# set_tests_properties(Anisotropic_Material_1D_Heat_Transfer
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
#                      
# add_test(NAME Anisotropic_Material_1D_Transient_Heat_Transfer
#          COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_anisotropic_transient_heat_transfer_material_1d)
# set_tests_properties(Anisotropic_Material_1D_Transient_Heat_Transfer
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
#          
#          
# # 2D Anisotropic Material Tests
# #-----------------------------------------------------------------------------
# add_test(NAME Anisotropic_Material_2D_Thermoelastic
#          COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_anisotropic_thermoelastic_material_2d)
# set_tests_properties(Anisotropic_Material_2D_Thermoelastic
#                     PROPERTIES 
#                     FIXTURES_REQUIRED   ConstantFieldFunction)
#                     
# add_test(NAME Anisotropic_Material_2D_Structural
#          COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_anisotropic_structural_material_2d)
# set_tests_properties(Anisotropic_Material_2D_Structural
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
# 
# add_test(NAME Anisotropic_Material_2D_Heat_Transfer
#          COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_anisotropic_heat_transfer_material_2d)
# set_tests_properties(Anisotropic_Material_2D_Heat_Transfer
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
#                     
# add_test(NAME Anisotropic_Material_2D_Transient_Heat_Transfer
#          COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_anisotropic_transient_heat_transfer_material_2d)
# set_tests_properties(Anisotropic_Material_2D_Transient_Heat_Transfer
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
#                     
#                     
# # 3D Anisotropic Material Tests
# #-----------------------------------------------------------------------------
# add_test(NAME Anisotropic_Material_3D_Thermoelastic
#          COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_anisotropic_thermoelastic_material_3d)
# set_tests_properties(Anisotropic_Material_3D_Thermoelastic
#                     PROPERTIES 
#                     FIXTURES_REQUIRED   ConstantFieldFunction)
# 
# add_test(NAME Anisotropic_Material_3D_Structural
#          COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_anisotropic_structural_material_3d)
# set_tests_properties(Anisotropic_Material_3D_Structural
#                      PROPERTIES 
#                      FIXTURES_REQUIRED  ConstantFieldFunction)
# 
# add_test(NAME Anisotropic_Material_3D_Heat_Transfer
#          COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_anisotropic_heat_transfer_material_3d)
# set_tests_properties(Anisotropic_Material_3D_Heat_Transfer
#                      PROPERTIES 
#                      FIXTURES_REQUIRED  ConstantFieldFunction)
#                      
# add_test(NAME Anisotropic_Material_3D_Transient_Heat_Transfer
#          COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_anisotropic_transient_heat_transfer_material_3d)
# set_tests_properties(Anisotropic_Material_3D_Transient_Heat_Transfer
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
# 
# add_test(NAME Anisotropic_Material_3D_Dynamic
#          COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests constant_anisotropic_dynamic_material_3d)
# set_tests_properties(Anisotropic_Material_3D_Dynamic
#                      PROPERTIES 
#                      FIXTURES_REQUIRED  ConstantFieldFunction)
#                      
