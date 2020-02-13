target_sources( mast_catch_tests
                PRIVATE
                ${CMAKE_CURRENT_LIST_DIR}/mast_isotropic_material.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_orthotropic_material.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_isotropic_material_sensitivity.cpp)

                
# ============================================================================
#                          ISOTROPIC MATERIAL TESTS
# ============================================================================
                
# 1D Isotropic Material Tests
#-----------------------------------------------------------------------------
add_test(NAME Isotropic_Material_1D_Thermoelastic
         COMMAND mast_catch_tests constant_isotropic_thermoelastic_material_1d)
set_tests_properties(Isotropic_Material_1D_Thermoelastic
                    PROPERTIES 
                    FIXTURES_REQUIRED   ConstantFieldFunction
                    FIXTURES_SETUP      Isotropic_Material_1D_Thermoelastic)

add_test(NAME Isotropic_Material_1D_Structural
         COMMAND mast_catch_tests constant_isotropic_structural_material_1d)
set_tests_properties(Isotropic_Material_1D_Structural
                    PROPERTIES 
                    FIXTURES_REQUIRED   ConstantFieldFunction
                    FIXTURES_SETUP      Isotropic_Material_1D_Structural)

add_test(NAME Isotropic_Material_1D_Heat_Transfer
         COMMAND mast_catch_tests constant_isotropic_heat_transfer_material_1d)
set_tests_properties(Isotropic_Material_1D_Heat_Transfer
                    PROPERTIES 
                    FIXTURES_REQUIRED   ConstantFieldFunction
                    FIXTURES_SETUP      Isotropic_Material_1D_Heat_Transfer)
                    
add_test(NAME Isotropic_Material_1D_Transient_Heat_Transfer
         COMMAND mast_catch_tests constant_isotropic_transient_heat_transfer_material_1d)
set_tests_properties(Isotropic_Material_1D_Transient_Heat_Transfer
                    PROPERTIES 
                    FIXTURES_REQUIRED   ConstantFieldFunction
                    FIXTURES_SETUP      Isotropic_Material_1D_Transient_Heat_Transfer)
                    
                    
# 2D Isotropic Material Tests
#-----------------------------------------------------------------------------
add_test(NAME Isotropic_Material_2D_Thermoelastic
         COMMAND mast_catch_tests constant_isotropic_thermoelastic_material_2d)
set_tests_properties(Isotropic_Material_2D_Thermoelastic
                    PROPERTIES 
                    FIXTURES_REQUIRED   ConstantFieldFunction
                    FIXTURES_SETUP      Isotropic_Material_2D_Thermoelastic)

add_test(NAME Isotropic_Material_2D_Structural
         COMMAND mast_catch_tests constant_isotropic_structural_material_2d)
set_tests_properties(Isotropic_Material_2D_Structural
                    PROPERTIES 
                    FIXTURES_REQUIRED   ConstantFieldFunction
                    FIXTURES_SETUP      Isotropic_Material_2D_Structural)

add_test(NAME Isotropic_Material_2D_Heat_Transfer
         COMMAND mast_catch_tests constant_isotropic_heat_transfer_material_2d)
set_tests_properties(Isotropic_Material_2D_Heat_Transfer
                     PROPERTIES 
                     FIXTURES_REQUIRED   ConstantFieldFunction
                     FIXTURES_SETUP      Isotropic_Material_2D_Heat_Transfer)
                    
add_test(NAME Isotropic_Material_2D_Transient_Heat_Transfer
         COMMAND mast_catch_tests constant_isotropic_transient_heat_transfer_material_2d)
set_tests_properties(Isotropic_Material_2D_Transient_Heat_Transfer
                     PROPERTIES 
                     FIXTURES_REQUIRED   ConstantFieldFunction
                     FIXTURES_SETUP      Isotropic_Material_2D_Transient_Heat_Transfer)
                    
                    
# 3D Isotropic Material Tests
#-----------------------------------------------------------------------------
add_test(NAME Isotropic_Material_3D_Thermoelastic
         COMMAND mast_catch_tests constant_isotropic_thermoelastic_material_3d)
set_tests_properties(Isotropic_Material_3D_Thermoelastic
                    PROPERTIES 
                    FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Isotropic_Material_3D_Structural
         COMMAND mast_catch_tests constant_isotropic_structural_material_3d)
set_tests_properties(Isotropic_Material_3D_Structural
                     PROPERTIES 
                     FIXTURES_REQUIRED  ConstantFieldFunction)

add_test(NAME Isotropic_Material_3D_Heat_Transfer
         COMMAND mast_catch_tests constant_isotropic_heat_transfer_material_3d)
set_tests_properties(Isotropic_Material_3D_Heat_Transfer
                     PROPERTIES 
                     FIXTURES_REQUIRED  ConstantFieldFunction)
                     
add_test(NAME Isotropic_Material_3D_Transient_Heat_Transfer
         COMMAND mast_catch_tests constant_isotropic_transient_heat_transfer_material_3d)
set_tests_properties(Isotropic_Material_3D_Transient_Heat_Transfer
                     PROPERTIES 
                     FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Isotropic_Material_3D_Dynamic
         COMMAND mast_catch_tests constant_isotropic_dynamic_material_3d)
set_tests_properties(Isotropic_Material_3D_Dynamic
                     PROPERTIES 
                     FIXTURES_REQUIRED  ConstantFieldFunction)
                    
                    
# ============================================================================
#                          ORTHOTROPIC MATERIAL TESTS
# ============================================================================

# 1D Orthotropic Material Tests
#-----------------------------------------------------------------------------
add_test(NAME Orthotropic_Material_1D_Thermoelastic
         COMMAND mast_catch_tests constant_orthotropic_thermoelastic_material_1d)
set_tests_properties(Orthotropic_Material_1D_Thermoelastic
                    PROPERTIES 
                    FIXTURES_REQUIRED   ConstantFieldFunction)
                    
add_test(NAME Orthotropic_Material_1D_Structural
         COMMAND mast_catch_tests constant_orthotropic_structural_material_1d)
set_tests_properties(Orthotropic_Material_1D_Structural
                     PROPERTIES 
                     FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Orthotropic_Material_1D_Heat_Transfer
         COMMAND mast_catch_tests constant_orthotropic_heat_transfer_material_1d)
set_tests_properties(Orthotropic_Material_1D_Heat_Transfer
                     PROPERTIES 
                     FIXTURES_REQUIRED   ConstantFieldFunction)
                    
add_test(NAME Orthotropic_Material_1D_Transient_Heat_Transfer
         COMMAND mast_catch_tests constant_orthotropic_transient_heat_transfer_material_1d)
set_tests_properties(Orthotropic_Material_1D_Transient_Heat_Transfer
                     PROPERTIES 
                     FIXTURES_REQUIRED   ConstantFieldFunction)

                     
# 2D Orthotropic Material Tests
#-----------------------------------------------------------------------------
add_test(NAME Orthotropic_Material_2D_Thermoelastic
         COMMAND mast_catch_tests constant_orthotropic_thermoelastic_material_2d)
set_tests_properties(Orthotropic_Material_2D_Thermoelastic
                    PROPERTIES 
                    FIXTURES_REQUIRED   ConstantFieldFunction)
                    
add_test(NAME Orthotropic_Material_2D_Structural
         COMMAND mast_catch_tests constant_orthotropic_structural_material_2d)
set_tests_properties(Orthotropic_Material_2D_Structural
                     PROPERTIES 
                     FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Orthotropic_Material_2D_Heat_Transfer
         COMMAND mast_catch_tests constant_orthotropic_heat_transfer_material_2d)
set_tests_properties(Orthotropic_Material_2D_Heat_Transfer
                     PROPERTIES 
                     FIXTURES_REQUIRED   ConstantFieldFunction)
                    
add_test(NAME Orthotropic_Material_2D_Transient_Heat_Transfer
         COMMAND mast_catch_tests constant_orthotropic_transient_heat_transfer_material_2d)
set_tests_properties(Orthotropic_Material_2D_Transient_Heat_Transfer
                     PROPERTIES 
                     FIXTURES_REQUIRED   ConstantFieldFunction)

                     
# 3D Orthotropic Material Tests
#-----------------------------------------------------------------------------
add_test(NAME Orthotropic_Material_3D_Thermoelastic
         COMMAND mast_catch_tests constant_orthotropic_thermoelastic_material_3d)
set_tests_properties(Orthotropic_Material_3D_Thermoelastic
                    PROPERTIES 
                    FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Orthotropic_Material_3D_Structural
         COMMAND mast_catch_tests constant_orthotropic_structural_material_3d)
set_tests_properties(Orthotropic_Material_3D_Structural
                     PROPERTIES 
                     FIXTURES_REQUIRED  ConstantFieldFunction)

add_test(NAME Orthotropic_Material_3D_Heat_Transfer
         COMMAND mast_catch_tests constant_orthotropic_heat_transfer_material_3d)
set_tests_properties(Orthotropic_Material_3D_Heat_Transfer
                     PROPERTIES 
                     FIXTURES_REQUIRED  ConstantFieldFunction)
                     
add_test(NAME Orthotropic_Material_3D_Transient_Heat_Transfer
         COMMAND mast_catch_tests constant_orthotropic_transient_heat_transfer_material_3d)
set_tests_properties(Orthotropic_Material_3D_Transient_Heat_Transfer
                     PROPERTIES 
                     FIXTURES_REQUIRED   ConstantFieldFunction)

add_test(NAME Orthotropic_Material_3D_Dynamic
         COMMAND mast_catch_tests constant_orthotropic_dynamic_material_3d)
set_tests_properties(Orthotropic_Material_3D_Dynamic
                     PROPERTIES 
                     FIXTURES_REQUIRED  ConstantFieldFunction)

         
# # ============================================================================
# #                          ANISOTROPIC MATERIAL TESTS
# # ============================================================================
# 
# # 1D Anisotropic Material Tests
# #-----------------------------------------------------------------------------
# add_test(NAME Anisotropic_Material_1D_Thermoelastic
#          COMMAND mast_catch_tests constant_anisotropic_thermoelastic_material_1d)
# set_tests_properties(Anisotropic_Material_1D_Thermoelastic
#                     PROPERTIES 
#                     FIXTURES_REQUIRED   ConstantFieldFunction)
#                     
# add_test(NAME Anisotropic_Material_1D_Structural
#          COMMAND mast_catch_tests constant_anisotropic_structural_material_1d)
# set_tests_properties(Anisotropic_Material_1D_Structural
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
# 
# add_test(NAME Anisotropic_Material_1D_Heat_Transfer
#          COMMAND mast_catch_tests constant_anisotropic_heat_transfer_material_1d)
# set_tests_properties(Anisotropic_Material_1D_Heat_Transfer
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
#                      
# add_test(NAME Anisotropic_Material_1D_Transient_Heat_Transfer
#          COMMAND mast_catch_tests constant_anisotropic_transient_heat_transfer_material_1d)
# set_tests_properties(Anisotropic_Material_1D_Transient_Heat_Transfer
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
#          
#          
# # 2D Anisotropic Material Tests
# #-----------------------------------------------------------------------------
# add_test(NAME Anisotropic_Material_2D_Thermoelastic
#          COMMAND mast_catch_tests constant_anisotropic_thermoelastic_material_2d)
# set_tests_properties(Anisotropic_Material_2D_Thermoelastic
#                     PROPERTIES 
#                     FIXTURES_REQUIRED   ConstantFieldFunction)
#                     
# add_test(NAME Anisotropic_Material_2D_Structural
#          COMMAND mast_catch_tests constant_anisotropic_structural_material_2d)
# set_tests_properties(Anisotropic_Material_2D_Structural
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
# 
# add_test(NAME Anisotropic_Material_2D_Heat_Transfer
#          COMMAND mast_catch_tests constant_anisotropic_heat_transfer_material_2d)
# set_tests_properties(Anisotropic_Material_2D_Heat_Transfer
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
#                     
# add_test(NAME Anisotropic_Material_2D_Transient_Heat_Transfer
#          COMMAND mast_catch_tests constant_anisotropic_transient_heat_transfer_material_2d)
# set_tests_properties(Anisotropic_Material_2D_Transient_Heat_Transfer
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
#                     
#                     
# # 3D Anisotropic Material Tests
# #-----------------------------------------------------------------------------
# add_test(NAME Anisotropic_Material_3D_Thermoelastic
#          COMMAND mast_catch_tests constant_anisotropic_thermoelastic_material_3d)
# set_tests_properties(Anisotropic_Material_3D_Thermoelastic
#                     PROPERTIES 
#                     FIXTURES_REQUIRED   ConstantFieldFunction)
# 
# add_test(NAME Anisotropic_Material_3D_Structural
#          COMMAND mast_catch_tests constant_anisotropic_structural_material_3d)
# set_tests_properties(Anisotropic_Material_3D_Structural
#                      PROPERTIES 
#                      FIXTURES_REQUIRED  ConstantFieldFunction)
# 
# add_test(NAME Anisotropic_Material_3D_Heat_Transfer
#          COMMAND mast_catch_tests constant_anisotropic_heat_transfer_material_3d)
# set_tests_properties(Anisotropic_Material_3D_Heat_Transfer
#                      PROPERTIES 
#                      FIXTURES_REQUIRED  ConstantFieldFunction)
#                      
# add_test(NAME Anisotropic_Material_3D_Transient_Heat_Transfer
#          COMMAND mast_catch_tests constant_anisotropic_transient_heat_transfer_material_3d)
# set_tests_properties(Anisotropic_Material_3D_Transient_Heat_Transfer
#                      PROPERTIES 
#                      FIXTURES_REQUIRED   ConstantFieldFunction)
# 
# add_test(NAME Anisotropic_Material_3D_Dynamic
#          COMMAND mast_catch_tests constant_anisotropic_dynamic_material_3d)
# set_tests_properties(Anisotropic_Material_3D_Dynamic
#                      PROPERTIES 
#                      FIXTURES_REQUIRED  ConstantFieldFunction)
#                      
