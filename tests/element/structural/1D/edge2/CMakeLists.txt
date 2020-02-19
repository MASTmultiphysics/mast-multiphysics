target_sources(mast_catch_tests
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_extension_internal_jacobian.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_extension_bending_internal_jacobian.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_extension_bending_coupling_internal_jacobian.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_extension_bending_shear_internal_jacobian.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_internal_jacobian.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_inertial_jacobian.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_thermal_jacobian.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_nonlinear_structural_thermal_jacobian.cpp)

        
# TODO: Implement patch test
# TODO: Implement surface pressure residual / Jacobian
# TODO: Implement piston theory residual / Jacobian testing
# TODO: Implement other elements (Edge3, etc.)
        
add_test(NAME Element_Edge2_Linear_Structural_Extension
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_extension_structural")
set_tests_properties(Element_Edge2_Linear_Structural_Extension
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED  Element_1D_Structural_Basic_Tests
        FIXTURES_SETUP  Element_Edge2_Linear_Structural_Extension)

add_test(NAME Element_Edge2_Linear_Structural_Extension_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_extension_structural")
set_tests_properties(Element_Edge2_Linear_Structural_Extension_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED  Element_1D_Structural_Basic_Tests_mpi
        FIXTURES_SETUP  Element_Edge2_Linear_Structural_Extension_mpi)



add_test(NAME Element_Edge2_Linear_Structural_Extension_Bending
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_extension_bending_structural")
set_tests_properties(Element_Edge2_Linear_Structural_Extension_Bending
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED  Element_Edge2_Linear_Structural_Extension
        FIXTURES_SETUP     Element_Edge2_Linear_Structural_Extension_Bending)

add_test(NAME Element_Edge2_Linear_Structural_Extension_Bending_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_extension_bending_structural")
set_tests_properties(Element_Edge2_Linear_Structural_Extension_Bending_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED  Element_Edge2_Linear_Structural_Extension_mpi
        FIXTURES_SETUP     Element_Edge2_Linear_Structural_Extension_Bending_mpi)



add_test(NAME Element_Edge2_Linear_Structural_Extension_Bending_Coupling
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_extension_bending_coupling_structural")
set_tests_properties(Element_Edge2_Linear_Structural_Extension_Bending_Coupling
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED  Element_Edge2_Linear_Structural_Extension_Bending
        FIXTURES_SETUP     Element_Edge2_Linear_Structural_Extension_Bending_Coupling)

add_test(NAME Element_Edge2_Linear_Structural_Extension_Bending_Coupling_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_extension_bending_coupling_structural")
set_tests_properties(Element_Edge2_Linear_Structural_Extension_Bending_Coupling_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED  Element_Edge2_Linear_Structural_Extension_Bending_mpi
        FIXTURES_SETUP     Element_Edge2_Linear_Structural_Extension_Bending_Coupling_mpi)



add_test(NAME Element_Edge2_Linear_Structural_Extension_Bending_Shear
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_extension_bending_shear_structural")
set_tests_properties(Element_Edge2_Linear_Structural_Extension_Bending_Shear
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED  Element_Edge2_Linear_Structural_Extension_Bending
        FIXTURES_SETUP     Element_Edge2_Linear_Structural_Extension_Bending_Shear)

add_test(NAME Element_Edge2_Linear_Structural_Extension_Bending_Shear_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_extension_bending_shear_structural")
set_tests_properties(Element_Edge2_Linear_Structural_Extension_Bending_Shear_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED  Element_Edge2_Linear_Structural_Extension_Bending_mpi
        FIXTURES_SETUP     Element_Edge2_Linear_Structural_Extension_Bending_Shear_mpi)



add_test(NAME Element_Edge2_Linear_Structural_All
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_structural")
set_tests_properties(Element_Edge2_Linear_Structural_All
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED "Element_Edge2_Linear_Structural_Extension_Bending_Coupling;Element_Edge2_Linear_Structural_Extension_Bending_Shear")

add_test(NAME Element_Edge2_Linear_Structural_All_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_structural")
set_tests_properties(Element_Edge2_Linear_Structural_All_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED "Element_Edge2_Linear_Structural_Extension_Bending_Coupling_mpi;Element_Edge2_Linear_Structural_Extension_Bending_Shear_mpi")



add_test(NAME Element_Edge2_Linear_Structural_Inertial_Consistent
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_structural_inertial_consistent")
set_tests_properties(Element_Edge2_Linear_Structural_Inertial_Consistent
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Element_1D_Structural_Basic_Tests)

add_test(NAME Element_Edge2_Linear_Structural_Inertial_Consistent_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_structural_inertial_consistent")
set_tests_properties(Element_Edge2_Linear_Structural_Inertial_Consistent_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Element_1D_Structural_Basic_Tests_mpi)



add_test(NAME Element_Edge2_Linear_Structural_Inertial_Lumped
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_structural_inertial_lumped")
set_tests_properties(Element_Edge2_Linear_Structural_Inertial_Lumped
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Element_1D_Structural_Basic_Tests)

add_test(NAME Element_Edge2_Linear_Structural_Inertial_Lumped_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_structural_inertial_lumped")
set_tests_properties(Element_Edge2_Linear_Structural_Inertial_Lumped_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Element_1D_Structural_Basic_Tests_mpi)



add_test(NAME Element_Edge2_Linear_Structural_Thermal
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_structural_thermal_jacobian")
set_tests_properties(Element_Edge2_Linear_Structural_Thermal
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED "Element_1D_Structural_Basic_Tests;Thermoelastic_Load")

add_test(NAME Element_Edge2_Linear_Structural_Thermal_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_linear_structural_thermal_jacobian")
set_tests_properties(Element_Edge2_Linear_Structural_Thermal_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED "Element_1D_Structural_Basic_Tests;Thermoelastic_Load_mpi")



add_test(NAME Element_Edge2_Nonlinear_Structural_Thermal
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_nonlinear_structural_thermal_jacobian")
set_tests_properties(Element_Edge2_Nonlinear_Structural_Thermal
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED "Element_1D_Structural_Basic_Tests;Thermoelastic_Load")

add_test(NAME Element_Edge2_Nonlinear_Structural_Thermal_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "edge2_nonlinear_structural_thermal_jacobian")
set_tests_properties(Element_Edge2_Nonlinear_Structural_Thermal_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED "Element_1D_Structural_Basic_Tests;Thermoelastic_Load_mpi")
