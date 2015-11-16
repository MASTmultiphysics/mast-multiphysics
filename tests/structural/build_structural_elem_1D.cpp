/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */



// MAST includes
#include "tests/structural/build_structural_elem_1D.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_element_base.h"


extern libMesh::LibMeshInit* _init;


MAST::BuildStructural1DElem::BuildStructural1DElem() {

    
    // create the mesh
    _mesh       = new libMesh::SerialMesh(_init->comm());
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_line(*_mesh, 1, 0, 2);
    
    // create the equation system
    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
    
    // create the libmesh system
    _sys       = &(_eq_sys->add_system<libMesh::NonlinearImplicitSystem>("structural"));
    
    // FEType to initialize the system
    libMesh::FEType fetype (libMesh::FIRST, libMesh::LAGRANGE);
    
    // initialize the system to the right set of variables
    _structural_sys = new MAST::StructuralSystemInitialization(*_sys,
                                                               _sys->name(),
                                                               fetype);
    _discipline     = new MAST::StructuralDiscipline(*_eq_sys);
    
    // initialize the equation system
    _eq_sys->init();
    
    // create the property functions and add them to the
    
    _thy             = new MAST::Parameter("thy",0.002);
    _thz             = new MAST::Parameter("thz",   1.);
    _E               = new MAST::Parameter("E",  72.e9);
    _nu              = new MAST::Parameter("nu",  0.33);
    _hy_off          = new MAST::Parameter("hyoff", 0.);
    _hz_off          = new MAST::Parameter("hzoff", 0.);
    _zero            = new MAST::Parameter("zero",  0.);
    
    

    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    _params_for_sensitivity.push_back(_E);
    _params_for_sensitivity.push_back(_nu);
    _params_for_sensitivity.push_back(_thy);
    _params_for_sensitivity.push_back(_thz);
    _params_for_sensitivity.push_back(_hy_off);
    _params_for_sensitivity.push_back(_hz_off);
    

    
    _thy_f           = new MAST::ConstantFieldFunction("hy",     *_thy);
    _thz_f           = new MAST::ConstantFieldFunction("hz",     *_thz);
    _E_f             = new MAST::ConstantFieldFunction("E",      *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",     *_nu);
    _hyoff_f         = new MAST::ConstantFieldFunction("hy_off", *_hy_off);
    _hzoff_f         = new MAST::ConstantFieldFunction("hz_off", *_hz_off);
    
    // create the material property card
    _m_card         = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(*_E_f);
    _m_card->add(*_nu_f);
    
    // create the element property card
    _p_card         = new MAST::Solid1DSectionElementPropertyCard;
    
    // tell the card about the orientation
    libMesh::Point orientation;
    orientation(1) = 1.;
    _p_card->y_vector() = orientation;
    
    // add the section properties to the card
    _p_card->add(*_thy_f);
    _p_card->add(*_thz_f);
    _p_card->add(*_hyoff_f);
    _p_card->add(*_hzoff_f);
    
    // tell the section property about the material property
    _p_card->set_material(*_m_card);
    
    _p_card->init();
    
    const unsigned int order = 1;
    Real
    mach    = 3.,
    a_inf   = 330.,
    gamma   = 1.4,
    rho     = 1.05;
    
    RealVectorX
    vel     = RealVectorX::Zero(3);
    
    // set velocity along x axis
    vel(0)  = 1.;
    
    // create the boundary condition
    _p_theory       = new MAST::PistonTheoryBoundaryCondition(order,
                                                              mach,
                                                              a_inf,
                                                              gamma,
                                                              rho,
                                                              vel);
    
    
}







MAST::BuildStructural1DElem::~BuildStructural1DElem() {
    
    delete _m_card;
    delete _p_card;
    
    delete _p_theory;
    
    delete _thy_f;
    delete _thz_f;
    delete _E_f;
    delete _nu_f;
    delete _hyoff_f;
    delete _hzoff_f;

    delete _thy;
    delete _thz;
    delete _hy_off;
    delete _hz_off;
    delete _E;
    delete _nu;
    delete _zero;
    
    
    

    delete _eq_sys;
    delete _mesh;
    
    delete _discipline;
    delete _structural_sys;
    
    
}


