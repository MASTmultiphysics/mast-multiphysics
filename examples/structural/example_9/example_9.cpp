// C/C++ Includes
#include <iostream>

// libMesh Includes
#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/point.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"

// MAST Includes
#include "base/nonlinear_system.h"
#include "elasticity/structural_system_initialization.h"
#include "base/physics_discipline_base.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/nonlinear_implicit_assembly.h"
#include "elasticity/structural_nonlinear_assembly.h"

#include "elasticity/structural_buckling_eigenproblem_assembly.h"
#include "elasticity/structural_buckling_eigenproblem_elem_operations.h"
#include "solver/slepc_eigen_solver.h"
#include <boundary_condition/point_load_condition.h>

#include "property_cards/solid_1d_bar_section_element_property_card.h"

#include "property_cards/cross_section_property_pilkey.h"

#include "libfort/fort.hpp"


// This is the main function which runs when the compiled code is called.
int main(int argc, const char** argv)
{
    // Initialize llibMesh Library
    libMesh::LibMeshInit init(argc, argv);
    
    const uint dim = 1;
        
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1420.5);         // Density
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    MAST::Parameter cp("cp_param",   908.0);          // Specific Heat Capacity
    MAST::Parameter k("k_param",     237.0);          // Thermal Conductivity
    
    const Real G = E() / (2. * (1. + nu()));
    
    // Define Section Properties as MAST Parameters
    MAST::Parameter DIM1("DIM1", 3.234);
    MAST::Parameter DIM2("DIM2", 1.422);
    MAST::Parameter offset_y("offy_param", 0.587);  // Section offset in y-direction
    MAST::Parameter offset_z("offz_param", -1.054); // Section offset in z-direction
    
    // Define Sensitivity Parameters
    std::vector<MAST::Parameter*> sens_params = {&DIM1, &DIM2};
    uint n_s = sens_params.size();
    
    // Create field functions to dsitribute these constant parameters throughout the model
    // Section Property Field Functions
    MAST::ConstantFieldFunction DIM1_f("DIM1", DIM1);
    MAST::ConstantFieldFunction DIM2_f("DIM2", DIM2);
    MAST::ConstantFieldFunction offsety_f("hy_off", offset_y);
    MAST::ConstantFieldFunction offsetz_f("hz_off", offset_z);
    // Material Property Field Functions
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    MAST::ConstantFieldFunction k_f("k_th", k);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(rho_f);
    material.add(k_f);
    material.add(cp_f);
    material.add(E_f);
    material.add(nu_f);
    
    // Initialize the section
    MAST::Solid1DBarSectionElementPropertyCard section;
    
    // Add the section property constant field functions to the section card
    section.add(DIM1_f);
    section.add(DIM2_f);
    section.add(offsety_f);
    section.add(offsetz_f);
    
    // Add the material card to the section card
    section.set_material(material);
    
    // Specify a section orientation point and add it to the section.
    RealVectorX orientation = RealVectorX::Zero(3);
    orientation(1) = 1.0;
    section.y_vector() = orientation;
    
    // Now initialize the section
    section.init(init);
    
    const libMesh::Point point(4.3, -3.5, -6.7);
    const Real time = 8.22;
    
    // Setup the table for section property console output
    fort::table properties_out;
    properties_out << fort::header << "Name" << "Property" << "Value" << fort::endr;
    
    fort::table dproperties_out;
    dproperties_out << fort::header << "Property" << "Value" << fort::endr;
    
    // Get Area
    const MAST::FieldFunction<Real>& Area = section.A();
    Real A, dA;
    Area(point, time, A);
    Area.derivative(DIM1, point, time, dA);
    properties_out << "Area" << "A" << std::to_string(A) << fort::endr;
    dproperties_out << "dA" << std::to_string(dA) << fort::endr;
    
    // Get Second Area Moments
    const MAST::FieldFunction<RealMatrixX>& SecondAreaMoments = section.I();
    Real Izz, Iyy, Izy;
    RealMatrixX I, dI;
    SecondAreaMoments(point, time, I);
    SecondAreaMoments.derivative(DIM1, point, time, dI);
    Izz = I(0,0);
    Iyy = I(1,1);
    Izy = I(0,1);
    properties_out << "Inertia" << "I_zz" << std::to_string(Izz) << fort::endr;
    properties_out << "Inertia" << "I_yy" << std::to_string(Iyy) << fort::endr;
    properties_out << "Inertia" << "I_zy" << std::to_string(Izy) << fort::endr;
    dproperties_out << "dI_zz" << std::to_string(dI(0,0)) << fort::endr;
    dproperties_out << "dI_yy" << std::to_string(dI(1,1)) << fort::endr;
    dproperties_out << "dI_zy" << std::to_string(dI(0,1)) << fort::endr;
    
    // Get Second Area Polar Moment
    const MAST::FieldFunction<Real>& PolarInteria = section.Ip();
    Real Ip, dIp;
    PolarInteria(point, time, Ip);
    PolarInteria.derivative(DIM1, point, time, dIp);
    properties_out << "Polar Inertia" << "I_xx" << std::to_string(Ip) << fort::endr;
    dproperties_out << "dI_xx" << std::to_string(dIp) << fort::endr;
    
    // Get Torsion Constant
    MAST::FieldFunction<Real>& TorsionConstant = section.J();
    Real J, dJ;
    TorsionConstant(point, time, J);
    TorsionConstant.derivative(DIM1, point, time, dJ);
    properties_out << "Torsion Constant" << "J" << std::to_string(J) << fort::endr;
    dproperties_out << "dJ" << std::to_string(dJ) << fort::endr;
    
    // Get first area moments
    const MAST::FieldFunction<Real>& AreaMomentY = section.Ay();
    const MAST::FieldFunction<Real>& AreaMomentZ = section.Az();
    Real Az, Ay, dAy, dAz;
    AreaMomentY(point, time, Ay);
    AreaMomentY.derivative(DIM1, point, time, dAy);
    AreaMomentZ(point, time, Az);
    AreaMomentZ.derivative(DIM1, point, time, dAz);
    properties_out << "First Area Moment" << "A_z" << std::to_string(Ay) << fort::endr;
    dproperties_out << "dA_z" << std::to_string(dAz) << fort::endr;
    properties_out << "First Area Moment" << "A_y" << std::to_string(Az) << fort::endr;
    dproperties_out << "dA_y" << std::to_string(dAy) << fort::endr;
    
    // Shear Coefficeints
    MAST::FieldFunction<RealMatrixX>& ShearCoefficients = section.Kap();
    RealMatrixX Kappa, dKappa;
    ShearCoefficients(point, time, Kappa);
    ShearCoefficients.derivative(DIM1, point, time, dKappa);
    properties_out << "Shear Coefficient" << "kappa_zz" << std::to_string(Kappa(0,0)) << fort::endr;
    properties_out << "Shear Coefficient" << "kappa_yy" << std::to_string(Kappa(1,1)) << fort::endr;
    properties_out << "Shear Coefficient" << "kappa_zy" << std::to_string(Kappa(0,1)) << fort::endr;
    dproperties_out << "dkappa_zz" << std::to_string(dKappa(0,0)) << fort::endr;
    dproperties_out << "dkappa_yy" << std::to_string(dKappa(1,1)) << fort::endr;
    dproperties_out << "dkappa_zy" << std::to_string(dKappa(0,1)) << fort::endr;
    
    // Get warping constant
    MAST::FieldFunction<Real>& WarpingConstant = section.Gam();
    Real W, dW;
    WarpingConstant(point, time, W);
    WarpingConstant.derivative(DIM1, point, time, dW);
    properties_out << "Warping Constant" << "W" << std::to_string(W) << fort::endr;
    dproperties_out << "dW" << std::to_string(dW) << fort::endr;
    
    libMesh::out << "\nProperty Values:\n" << properties_out.to_string() << std::endl;
    
    libMesh::out << "Property Derivative w.r.t. DIM1\n" << dproperties_out.to_string() << std::endl;
    
    return 0;
}
