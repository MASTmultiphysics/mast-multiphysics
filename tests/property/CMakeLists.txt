target_sources(mast_catch_tests
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/mast_solid_1d_section_element_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_solid_1d_arbitrary_section_element_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_solid_1d_rod_section_element_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_solid_2d_section_element_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_3d_isotropic_element_property_card.cpp)
       

## ============================================================================
##                              3D Property Card Tests
## ============================================================================
# Thermoelastic
add_test(NAME Element_Property_Card_Isotropic_3D_Thermoelastic
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_thermoelastic_isotropic_3d")
set_tests_properties(Element_Property_Card_Isotropic_3D_Thermoelastic
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Isotropic_Material_3D_Thermoelastic
        FIXTURES_SETUP Element_Property_Card_Isotropic_3D_Thermoelastic)

add_test(NAME Element_Property_Card_Isotropic_3D_Thermoelastic_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_thermoelastic_isotropic_3d")
set_tests_properties(Element_Property_Card_Isotropic_3D_Thermoelastic_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Isotropic_Material_3D_Thermoelastic_mpi
        FIXTURES_SETUP Element_Property_Card_Isotropic_3D_Thermoelastic_mpi)


# Structural
add_test(NAME Element_Property_Card_Isotropic_3D_Structural
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_structural_isotropic_3d")
set_tests_properties(Element_Property_Card_Isotropic_3D_Structural
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Isotropic_Material_3D_Structural
        FIXTURES_SETUP Element_Property_Card_Isotropic_3D_Structural)

add_test(NAME Element_Property_Card_Isotropic_3D_Structural_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_structural_isotropic_3d")
set_tests_properties(Element_Property_Card_Isotropic_3D_Structural_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Isotropic_Material_3D_Structural_mpi
        FIXTURES_SETUP Element_Property_Card_Isotropic_3D_Structural_mpi)


# Dynamic
add_test(NAME Element_Property_Card_Isotropic_3D_Dynamic
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_dynamic_isotropic_3d")
set_tests_properties(Element_Property_Card_Isotropic_3D_Dynamic
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Isotropic_Material_3D_Dynamic
        FIXTURES_SETUP Element_Property_Card_Isotropic_3D_Dynamic)

add_test(NAME Element_Property_Card_Isotropic_3D_Dynamic_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_dynamic_isotropic_3d")
set_tests_properties(Element_Property_Card_Isotropic_3D_Dynamic_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Isotropic_Material_3D_Dynamic_mpi
        FIXTURES_SETUP Element_Property_Card_Isotropic_3D_Dynamic_mpi)


# Heat Transfer
add_test(NAME Element_Property_Card_Isotropic_3D_Heat_Transfer
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_heat_transfer_isotropic_3d")
set_tests_properties(Element_Property_Card_Isotropic_3D_Heat_Transfer
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Isotropic_Material_3D_Heat_Transfer
        FIXTURES_SETUP Element_Property_Card_3D_Heat_Transfer)

add_test(NAME Element_Property_Card_Isotropic_3D_Heat_Transfer_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_heat_transfer_isotropic_3d")
set_tests_properties(Element_Property_Card_Isotropic_3D_Heat_Transfer_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Isotropic_Material_3D_Heat_Transfer_mpi
        FIXTURES_SETUP Element_Property_Card_3D_Heat_Transfer_mpi)


# Transient Heat Transfer
add_test(NAME Element_Property_Card_Isotropic_3D_Transient_Heat_Transfer
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_transient_heat_transfer_isotropic_3d")
set_tests_properties(Element_Property_Card_Isotropic_3D_Transient_Heat_Transfer
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Isotropic_Material_3D_Transient_Heat_Transfer
        FIXTURES_SETUP Element_Property_Card_3D_Transient_Heat_Transfer)

add_test(NAME Element_Property_Card_Isotropic_3D_Transient_Heat_Transfer_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_transient_heat_transfer_isotropic_3d")
set_tests_properties(Element_Property_Card_Isotropic_3D_Transient_Heat_Transfer_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Isotropic_Material_3D_Transient_Heat_Transfer_mpi
        FIXTURES_SETUP Element_Property_Card_3D_Transient_Heat_Transfer_mpi)

## ============================================================================
##                              2D Property Card Tests
## ============================================================================
# Structural
add_test(NAME Element_Property_Card_2D_Structural
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_structural_2d")
set_tests_properties(Element_Property_Card_2D_Structural
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Isotropic_Material_2D_Structural
        FIXTURES_SETUP Element_Property_Card_2D_Structural)

add_test(NAME Element_Property_Card_2D_Structural_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_structural_2d")
set_tests_properties(Element_Property_Card_2D_Structural_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Isotropic_Material_2D_Structural_mpi
        FIXTURES_SETUP Element_Property_Card_2D_Structural_mpi)


# Thermoelastic
add_test(NAME Element_Property_Card_2D_Thermoelastic
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_thermoelastic_2d")
set_tests_properties(Element_Property_Card_2D_Thermoelastic
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Isotropic_Material_2D_Thermoelastic
        FIXTURES_SETUP Element_Property_Card_2D_Thermoelastic)

add_test(NAME Element_Property_Card_2D_Thermoelastic_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_thermoelastic_2d")
set_tests_properties(Element_Property_Card_2D_Thermoelastic_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Isotropic_Material_2D_Thermoelastic_mpi
        FIXTURES_SETUP Element_Property_Card_2D_Thermoelastic_mpi)


# Dynamic
add_test(NAME Element_Property_Card_2D_Dynamic
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_dynamic_2d")
set_tests_properties(Element_Property_Card_2D_Dynamic
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_SETUP Element_Property_Card_2D_Dynamic)

add_test(NAME Element_Property_Card_2D_Dynamic_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_dynamic_2d")
set_tests_properties(Element_Property_Card_2D_Dynamic_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_SETUP Element_Property_Card_2D_Dynamic_mpi)


# Heat Transfer
add_test(NAME Element_Property_Card_2D_Heat_Transfer
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_heat_transfer_2d")
set_tests_properties(Element_Property_Card_2D_Heat_Transfer
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Isotropic_Material_2D_Heat_Transfer
        FIXTURES_SETUP Element_Property_Card_2D_Heat_Transfer)

add_test(NAME Element_Property_Card_2D_Heat_Transfer_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_heat_transfer_2d")
set_tests_properties(Element_Property_Card_2D_Heat_Transfer_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Isotropic_Material_2D_Heat_Transfer_mpi
        FIXTURES_SETUP Element_Property_Card_2D_Heat_Transfer_mpi)


# Transient Heat Transfer
add_test(NAME Element_Property_Card_2D_Transient_Heat_Transfer
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_transient_heat_transfer_2d")
set_tests_properties(Element_Property_Card_2D_Transient_Heat_Transfer
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Isotropic_Material_2D_Transient_Heat_Transfer
        FIXTURES_SETUP Element_Property_Card_2D_Transient_Heat_Transfer)

add_test(NAME Element_Property_Card_2D_Transient_Heat_Transfer_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "element_property_card_constant_transient_heat_transfer_2d")
set_tests_properties(Element_Property_Card_2D_Transient_Heat_Transfer_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Isotropic_Material_2D_Transient_Heat_Transfer_mpi
        FIXTURES_SETUP Element_Property_Card_2D_Transient_Heat_Transfer_mpi)


## ============================================================================
##                              1D Property Card Tests
## ============================================================================
     
## ============================================================================
# 1D (Default, Rectangular) Property Card Tests
add_test(NAME Element_Property_Card_1D_Base
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "solid_element_property_card_constant_base_1d")
set_tests_properties(Element_Property_Card_1D_Base
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Isotropic_Material_1D_Structural)

add_test(NAME Element_Property_Card_1D_Base_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "solid_element_property_card_constant_base_1d")
set_tests_properties(Element_Property_Card_1D_Base_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Isotropic_Material_1D_Structural_mpi)

                     
add_test(NAME Element_Property_Card_1D_Base_Sensitivity
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "solid_element_property_card_constant_base_sensitivity_1d")
set_tests_properties(Element_Property_Card_1D_Base_Sensitivity
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Element_Property_Card_1D_Base)

add_test(NAME Element_Property_Card_1D_Base_Sensitivity_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "solid_element_property_card_constant_base_sensitivity_1d")
set_tests_properties(Element_Property_Card_1D_Base_Sensitivity_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Element_Property_Card_1D_Base_mpi)


add_test(NAME Element_Property_Card_1D_Heat_Transfer
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "solid_element_property_card_constant_heat_transfer_1d")
set_tests_properties(Element_Property_Card_1D_Heat_Transfer
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Isotropic_Material_1D_Heat_Transfer)

add_test(NAME Element_Property_Card_1D_Heat_Transfer_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "solid_element_property_card_constant_heat_transfer_1d")
set_tests_properties(Element_Property_Card_1D_Heat_Transfer_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Isotropic_Material_1D_Heat_Transfer_mpi)


add_test(NAME Element_Property_Card_1D_Thermoelastic
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "solid_element_property_card_constant_thermoelastic_1d")
set_tests_properties(Element_Property_Card_1D_Thermoelastic
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Isotropic_Material_1D_Thermoelastic)

add_test(NAME Element_Property_Card_1D_Thermoelastic_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "solid_element_property_card_constant_thermoelastic_1d")
set_tests_properties(Element_Property_Card_1D_Thermoelastic_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Isotropic_Material_1D_Thermoelastic_mpi)


add_test(NAME Element_Property_Card_1D_Structural
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "solid_element_property_card_constant_structural_1d")
set_tests_properties(Element_Property_Card_1D_Structural
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED  Isotropic_Material_1D_Structural)

add_test(NAME Element_Property_Card_1D_Structural_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "solid_element_property_card_constant_structural_1d")
set_tests_properties(Element_Property_Card_1D_Structural_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED  Isotropic_Material_1D_Structural_mpi)


add_test(NAME Element_Property_Card_1D_Dynamic
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "solid_element_property_card_constant_dynamic_1d")
set_tests_properties(Element_Property_Card_1D_Dynamic
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Isotropic_Material_1D_Dynamic)

add_test(NAME Element_Property_Card_1D_Dynamic_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "solid_element_property_card_constant_dynamic_1d")
set_tests_properties(Element_Property_Card_1D_Dynamic_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Isotropic_Material_1D_Dynamic_mpi)


## ============================================================================
# 1D Arbitrary Cross Section
add_test(NAME Arbitrary_Element_Property_Card_1D_Base
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "arbitrary_element_property_card_constant_base_1d")
set_tests_properties(Arbitrary_Element_Property_Card_1D_Base
                     PROPERTIES
                     LABELS "SEQ"
                     FIXTURES_REQUIRED Isotropic_Material_1D_Structural)
                     
add_test(NAME Arbitrary_Element_Property_Card_1D_Base_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "arbitrary_element_property_card_constant_base_1d")
set_tests_properties(Arbitrary_Element_Property_Card_1D_Base_mpi
                     PROPERTIES
                     LABELS "MPI"
                     FIXTURES_REQUIRED Isotropic_Material_1D_Structural_mpi)
## ============================================================================


# ## ============================================================================
# # 1D BAR Cross Section
# add_test(NAME Bar_Element_Property_Card_1D_Base
#          COMMAND mast_catch_tests "bar_element_property_card_constant_base_1d")
# set_tests_properties(Bar_Element_Property_Card_1D_Base
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Structural)
#                      
# add_test(NAME Bar_Element_Property_Card_1D_Base_Sensitivity
#          COMMAND mast_catch_tests "bar_element_property_card_constant_base_sensitivity_1d")
# set_tests_properties(Bar_Element_Property_Card_1D_Base_Sensitivity
#                      PROPERTIES
#                      FIXTURES_REQUIRED Bar_Element_Property_Card_1D_Base)
#                      
# add_test(NAME Bar_Element_Property_Card_1D_Heat_Transfer
#          COMMAND mast_catch_tests "bar_element_property_card_constant_heat_transfer_1d")
# set_tests_properties(Bar_Element_Property_Card_1D_Heat_Transfer
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Heat_Transfer)
# 
# add_test(NAME Bar_Element_Property_Card_1D_Thermoelastic
#          COMMAND mast_catch_tests "bar_element_property_card_constant_thermoelastic_1d")
# set_tests_properties(Bar_Element_Property_Card_1D_Thermoelastic
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Thermoelastic)
# 
# add_test(NAME Bar_Element_Property_Card_1D_Structural
#          COMMAND mast_catch_tests "bar_element_property_card_constant_structural_1d")
# set_tests_properties(Bar_Element_Property_Card_1D_Structural
#                      PROPERTIES
#                      FIXTURES_REQUIRED  Isotropic_Material_1D_Structural)
#                      
# add_test(NAME Bar_Element_Property_Card_1D_Dynamic
#          COMMAND mast_catch_tests "bar_element_property_card_constant_dynamic_1d")
# set_tests_properties(Bar_Element_Property_Card_1D_Dynamic
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Dynamic)
# ## ============================================================================


## ============================================================================
# 1D ROD Cross Section
add_test(NAME Rod_Element_Property_Card_1D_Base
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "rod_element_property_card_constant_base_1d")
set_tests_properties(Rod_Element_Property_Card_1D_Base
                     PROPERTIES
                     LABELS "SEQ"
                     FIXTURES_REQUIRED Isotropic_Material_1D_Structural)
                     
add_test(NAME Rod_Element_Property_Card_1D_Base_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "rod_element_property_card_constant_base_1d")
set_tests_properties(Rod_Element_Property_Card_1D_Base_mpi
                     PROPERTIES
                     LABELS "MPI"
                     FIXTURES_REQUIRED Isotropic_Material_1D_Structural_mpi)

add_test(NAME Rod_Element_Property_Card_1D_Base_Sensitivity
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "rod_element_property_card_constant_base_sensitivity_1d")
set_tests_properties(Rod_Element_Property_Card_1D_Base_Sensitivity
                     PROPERTIES
                     LABELS "SEQ"
                     FIXTURES_REQUIRED Rod_Element_Property_Card_1D_Base)
                     
add_test(NAME Rod_Element_Property_Card_1D_Heat_Transfer
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "rod_element_property_card_constant_heat_transfer_1d")
set_tests_properties(Rod_Element_Property_Card_1D_Heat_Transfer
                     PROPERTIES
                     LABELS "SEQ"
                     FIXTURES_REQUIRED Isotropic_Material_1D_Heat_Transfer)

add_test(NAME Rod_Element_Property_Card_1D_Thermoelastic
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "rod_element_property_card_constant_thermoelastic_1d")
set_tests_properties(Rod_Element_Property_Card_1D_Thermoelastic
                     PROPERTIES
                     LABELS "SEQ"
                     FIXTURES_REQUIRED Isotropic_Material_1D_Thermoelastic)

add_test(NAME Rod_Element_Property_Card_1D_Structural
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "rod_element_property_card_constant_structural_1d")
set_tests_properties(Rod_Element_Property_Card_1D_Structural
                     PROPERTIES
                     LABELS "SEQ"
                     FIXTURES_REQUIRED  Isotropic_Material_1D_Structural)
                     
add_test(NAME Rod_Element_Property_Card_1D_Dynamic
         COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "rod_element_property_card_constant_dynamic_1d")
set_tests_properties(Rod_Element_Property_Card_1D_Dynamic
                     PROPERTIES
                     LABELS "SEQ"
                     FIXTURES_REQUIRED Isotropic_Material_1D_Dynamic)
## ============================================================================


# ## ============================================================================
# # 1D TUBE Cross Section
# add_test(NAME Tube_Element_Property_Card_1D_Base
#          COMMAND mast_catch_tests "tube_element_property_card_constant_base_1d")
# set_tests_properties(Tube_Element_Property_Card_1D_Base
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Structural)
#                      
# add_test(NAME Tube_Element_Property_Card_1D_Heat_Transfer
#          COMMAND mast_catch_tests "tube_element_property_card_constant_heat_transfer_1d")
# set_tests_properties(Tube_Element_Property_Card_1D_Heat_Transfer
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Heat_Transfer)
# 
# add_test(NAME Tube_Element_Property_Card_1D_Thermoelastic
#          COMMAND mast_catch_tests "tube_element_property_card_constant_thermoelastic_1d")
# set_tests_properties(Tube_Element_Property_Card_1D_Thermoelastic
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Thermoelastic)
# 
# add_test(NAME Tube_Element_Property_Card_1D_Structural
#          COMMAND mast_catch_tests "tube_element_property_card_constant_structural_1d")
# set_tests_properties(Tube_Element_Property_Card_1D_Structural
#                      PROPERTIES
#                      FIXTURES_REQUIRED  Isotropic_Material_1D_Structural)
#                      
# add_test(NAME Tube_Element_Property_Card_1D_Dynamic
#          COMMAND mast_catch_tests "tube_element_property_card_constant_dynamic_1d")
# set_tests_properties(Tube_Element_Property_Card_1D_Dynamic
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Dynamic)
# ## ============================================================================


# ## ============================================================================
# # 1D TUBE2 Cross Section
# add_test(NAME Tube2_Element_Property_Card_1D_Base
#          COMMAND mast_catch_tests "tube2_element_property_card_constant_base_1d")
# set_tests_properties(Tube2_Element_Property_Card_1D_Base
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Structural)
#                      
# add_test(NAME Tube2_Element_Property_Card_1D_Heat_Transfer
#          COMMAND mast_catch_tests "tube2_element_property_card_constant_heat_transfer_1d")
# set_tests_properties(Tube2_Element_Property_Card_1D_Heat_Transfer
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Heat_Transfer)
# 
# add_test(NAME Tube2_Element_Property_Card_1D_Thermoelastic
#          COMMAND mast_catch_tests "tube2_element_property_card_constant_thermoelastic_1d")
# set_tests_properties(Tube2_Element_Property_Card_1D_Thermoelastic
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Thermoelastic)
# 
# add_test(NAME Tube2_Element_Property_Card_1D_Structural
#          COMMAND mast_catch_tests "tube2_element_property_card_constant_structural_1d")
# set_tests_properties(Tube2_Element_Property_Card_1D_Structural
#                      PROPERTIES
#                      FIXTURES_REQUIRED  Isotropic_Material_1D_Structural)
#                      
# add_test(NAME Tube2_Element_Property_Card_1D_Dynamic
#          COMMAND mast_catch_tests "tube2_element_property_card_constant_dynamic_1d")
# set_tests_properties(Tube2_Element_Property_Card_1D_Dynamic
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Dynamic)
# ## ============================================================================


# ## ============================================================================
# # 1D I1 Cross Section
# add_test(NAME I1_Element_Property_Card_1D_Base
#          COMMAND mast_catch_tests "I1_element_property_card_constant_base_1d")
# set_tests_properties(I1_Element_Property_Card_1D_Base
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Structural)
#                      
# add_test(NAME I1_Element_Property_Card_1D_Base_Sensitivity
#          COMMAND mast_catch_tests "I1_element_property_card_constant_base_sensitivity_1d")
# set_tests_properties(I1_Element_Property_Card_1D_Base_Sensitivity
#                      PROPERTIES
#                      FIXTURES_REQUIRED I1_Element_Property_Card_1D_Base)
#                      
# add_test(NAME I1_Element_Property_Card_1D_Heat_Transfer
#          COMMAND mast_catch_tests "I1_element_property_card_constant_heat_transfer_1d")
# set_tests_properties(I1_Element_Property_Card_1D_Heat_Transfer
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Heat_Transfer)
# 
# add_test(NAME I1_Element_Property_Card_1D_Thermoelastic
#          COMMAND mast_catch_tests "I1_element_property_card_constant_thermoelastic_1d")
# set_tests_properties(I1_Element_Property_Card_1D_Thermoelastic
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Thermoelastic)
# 
# add_test(NAME I1_Element_Property_Card_1D_Structural
#          COMMAND mast_catch_tests "I1_element_property_card_constant_structural_1d")
# set_tests_properties(I1_Element_Property_Card_1D_Structural
#                      PROPERTIES
#                      FIXTURES_REQUIRED  Isotropic_Material_1D_Structural)
#                      
# add_test(NAME I1_Element_Property_Card_1D_Dynamic
#          COMMAND mast_catch_tests "I1_element_property_card_constant_dynamic_1d")
# ## ============================================================================
                     

# ## ============================================================================
# # 1D L Cross Section
# add_test(NAME L_Element_Property_Card_1D_Base
#          COMMAND mast_catch_tests "L_element_property_card_constant_base_1d")
# set_tests_properties(L_Element_Property_Card_1D_Base
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Structural)
#                      
# add_test(NAME L_Element_Property_Card_1D_Heat_Transfer
#          COMMAND mast_catch_tests "L_element_property_card_constant_heat_transfer_1d")
# set_tests_properties(L_Element_Property_Card_1D_Heat_Transfer
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Heat_Transfer)
# 
# add_test(NAME L_Element_Property_Card_1D_Thermoelastic
#          COMMAND mast_catch_tests "L_element_property_card_constant_thermoelastic_1d")
# set_tests_properties(L_Element_Property_Card_1D_Thermoelastic
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Thermoelastic)
# 
# add_test(NAME L_Element_Property_Card_1D_Structural
#          COMMAND mast_catch_tests "L_element_property_card_constant_structural_1d")
# set_tests_properties(L_Element_Property_Card_1D_Structural
#                      PROPERTIES
#                      FIXTURES_REQUIRED  Isotropic_Material_1D_Structural)
#                      
# add_test(NAME L_Element_Property_Card_1D_Dynamic
#          COMMAND mast_catch_tests "L_element_property_card_constant_dynamic_1d")
# set_tests_properties(L_Element_Property_Card_1D_Dynamic
#                      PROPERTIES
#                      FIXTURES_REQUIRED Isotropic_Material_1D_Dynamic)
# ## ============================================================================
