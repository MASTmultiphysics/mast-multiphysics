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
#include "tests/structural/build_structural_elem_2D.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_element_base.h"


extern libMesh::LibMeshInit* _init;


MAST::BuildStructural2DElem::BuildStructural2DElem():
_initialized(false),
_e_type(libMesh::INVALID_ELEM),
_mesh(NULL),
_eq_sys(NULL),
_sys(NULL),
_structural_sys(NULL),
_discipline(NULL),
_thz(NULL),
_E(NULL),
_nu(NULL),
_kappa(NULL),
_hzoff(NULL),
_zero(NULL),
_thz_f(NULL),
_E_f(NULL),
_nu_f(NULL),
_kappa_f(NULL),
_hzoff_f(NULL),
_m_card(NULL),
_p_card(NULL),
_p_theory(NULL) {
    
}



void
MAST::BuildStructural2DElem::init(bool if_link_offset_to_th,
                                  bool if_nonlinear,
                                  libMesh::ElemType e_type) {
    
    // make sure that this has not already been initialized
    libmesh_assert(!_initialized);
    _e_type = e_type;
    
    // create the mesh
    _mesh       = new libMesh::SerialMesh(_init->comm());
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_square(*_mesh,
                                                 1, 1,
                                                 0, 2,
                                                 0, 2,
                                                 e_type);
    
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
    
    _thz             = new MAST::Parameter(  "thz",    0.002);
    _E               = new MAST::Parameter(    "E",    72.e9);
    _nu              = new MAST::Parameter(   "nu",     0.33);
    _kappa           = new MAST::Parameter("kappa",    5./6.);
    _zero            = new MAST::Parameter( "zero",       0.);
    _hzoff           = new MAST::Parameter( "hzoff",      0.);
    _temp            = new MAST::Parameter("temp",       60.);
    _alpha           = new MAST::Parameter("alpha",   2.5e-5);
    
    
    
    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    _params_for_sensitivity.push_back(    _E);
    _params_for_sensitivity.push_back(   _nu);
    _params_for_sensitivity.push_back(  _thz);
    _params_for_sensitivity.push_back(_kappa);
    
    
    
    _thz_f           = new MAST::ConstantFieldFunction(    "h",     *_thz);
    _E_f             = new MAST::ConstantFieldFunction(    "E",       *_E);
    _nu_f            = new MAST::ConstantFieldFunction(   "nu",      *_nu);
    _kappa_f         = new MAST::ConstantFieldFunction("kappa",   *_kappa);
    _temp_f          = new MAST::ConstantFieldFunction("temperature", *_temp);
    _ref_temp_f      = new MAST::ConstantFieldFunction("ref_temperature", *_zero);
    _alpha_f         = new MAST::ConstantFieldFunction("alpha_expansion", *_alpha);
    if (!if_link_offset_to_th)
        _hzoff_f         = new MAST::ConstantFieldFunction(  "off",   *_hzoff);
    else
        _hzoff_f         = new MAST::ConstantFieldFunction(  "off",   *_thz);
    
    // create the material property card
    _m_card         = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(    *_E_f);
    _m_card->add(   *_nu_f);
    _m_card->add(*_kappa_f);
    _m_card->add(*_alpha_f);
    
    // create the element property card
    _p_card         = new MAST::Solid2DSectionElementPropertyCard;
    
    // add the section properties to the card
    _p_card->add(*_thz_f);
    _p_card->add(*_hzoff_f);
    
    // tell the section property about the material property
    _p_card->set_material(*_m_card);
    if (if_nonlinear) _p_card->set_strain(MAST::VON_KARMAN_STRAIN);
    
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
    
    
    _thermal_load   = new MAST::BoundaryConditionBase(MAST::TEMPERATURE);
    _thermal_load->add(*_temp_f);
    _thermal_load->add(*_ref_temp_f);
}







MAST::BuildStructural2DElem::~BuildStructural2DElem() {
    
    delete _m_card;
    delete _p_card;
    
    delete _p_theory;
    delete _thermal_load;
    
    delete _thz_f;
    delete _E_f;
    delete _nu_f;
    delete _hzoff_f;
    delete _temp_f;
    delete _ref_temp_f;
    delete _alpha_f;
    
    delete _thz;
    delete _hzoff;
    delete _E;
    delete _nu;
    delete _zero;
    delete _temp;
    delete _alpha;
    
    
    
    delete _eq_sys;
    delete _mesh;
    
    delete _discipline;
    delete _structural_sys;
    
    
}


