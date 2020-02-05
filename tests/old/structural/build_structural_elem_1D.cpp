/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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
#include "base/nonlinear_system.h"

extern libMesh::LibMeshInit* __init;


MAST::BuildStructural1DElem::BuildStructural1DElem():
_initialized(false),
_mesh(nullptr),
_eq_sys(nullptr),
_sys(nullptr),
_structural_sys(nullptr),
_discipline(nullptr),
_thy(nullptr),
_thz(nullptr),
_E(nullptr),
_nu(nullptr),
_hy_off(nullptr),
_hz_off(nullptr),
_zero(nullptr),
_temp(nullptr),
_alpha(nullptr),
_thy_f(nullptr),
_thz_f(nullptr),
_E_f(nullptr),
_nu_f(nullptr),
_hyoff_f(nullptr),
_hzoff_f(nullptr),
_temp_f(nullptr),
_ref_temp_f(nullptr),
_alpha_f(nullptr),
_m_card(nullptr),
_p_card(nullptr),
_p_theory(nullptr),
_thermal_load(nullptr) {
    
 }


void
MAST::BuildStructural1DElem::init(bool if_link_offset_to_th,
                                  bool if_nonlinear) {

    // make sure that this has not already been initialized
    libmesh_assert(!_initialized);

    // create the mesh
    _mesh       = new libMesh::SerialMesh(__init->comm());
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_line(*_mesh, 1, 0, 2);
    
    // create the equation system
    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
    
    // create the libmesh system
    _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
    
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
    
    _thy             = new MAST::Parameter("thy",      0.06);
    _thz             = new MAST::Parameter("thz",      0.02);
    _E               = new MAST::Parameter("E",       72.e9);
    _nu              = new MAST::Parameter("nu",       0.33);
    _hy_off          = new MAST::Parameter("hyoff",      0.);
    _hz_off          = new MAST::Parameter("hzoff",      0.);
    _zero            = new MAST::Parameter("zero",       0.);
    _temp            = new MAST::Parameter("temp",      60.);
    _alpha           = new MAST::Parameter("alpha",  2.5e-5);
    _velocity        = new MAST::Parameter("V"   ,      200.);
    _mach            = new MAST::Parameter("mach",        4.);
    _rho_air         = new MAST::Parameter("rho" ,       1.1);
    _gamma_air       = new MAST::Parameter("gamma",      1.4);
    _dwdx            = new MAST::Parameter("dwdx",       2.5);
    _dwdt            = new MAST::Parameter("dwdx",       4.5);
    
    

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
    _temp_f          = new MAST::ConstantFieldFunction("temperature", *_temp);
    _ref_temp_f      = new MAST::ConstantFieldFunction("ref_temperature", *_zero);
    _alpha_f         = new MAST::ConstantFieldFunction("alpha_expansion", *_alpha);
    if (!if_link_offset_to_th) {
        _hyoff_f         = new MAST::ConstantFieldFunction("hy_off", *_hy_off);
        _hzoff_f         = new MAST::ConstantFieldFunction("hz_off", *_hz_off);
    }
    else {
        _hyoff_f         = new MAST::ConstantFieldFunction("hy_off", *_thy);
        _hzoff_f         = new MAST::ConstantFieldFunction("hz_off", *_thz);
    }
    _velocity_f      = new MAST::ConstantFieldFunction("V",      *_velocity);
    _mach_f          = new MAST::ConstantFieldFunction("mach",       *_mach);
    _rho_air_f       = new MAST::ConstantFieldFunction("rho",     *_rho_air);
    _gamma_air_f     = new MAST::ConstantFieldFunction("gamma", *_gamma_air);
    _dwdx_f          = new MAST::ConstantFieldFunction("dwdx",       *_dwdx);
    _dwdt_f          = new MAST::ConstantFieldFunction("dwdt",      *_dwdt);

    
    // create the material property card
    _m_card         = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(*_E_f);
    _m_card->add(*_nu_f);
    _m_card->add(*_alpha_f);
    
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
    if (if_nonlinear) _p_card->set_strain(MAST::VON_KARMAN_STRAIN);
    
    _p_card->init();
    
    RealVectorX
    vel     = RealVectorX::Zero(3);
    
    // set velocity along x axis
    vel(0)  = 1.;
    
    // create the boundary condition
    _p_theory       = new MAST::PistonTheoryBoundaryCondition(3,     // order
                                                              vel);
    _p_theory->add(*_velocity_f);
    _p_theory->add(*_mach_f);
    _p_theory->add(*_rho_air_f);
    _p_theory->add(*_gamma_air_f);
    
    _thermal_load   = new MAST::BoundaryConditionBase(MAST::TEMPERATURE);
    _thermal_load->add(*_temp_f);
    _thermal_load->add(*_ref_temp_f);
}







MAST::BuildStructural1DElem::~BuildStructural1DElem() {
    
    delete _m_card;
    delete _p_card;
    
    delete _p_theory;
    delete _thermal_load;
    
    delete _thy_f;
    delete _thz_f;
    delete _E_f;
    delete _nu_f;
    delete _hyoff_f;
    delete _hzoff_f;
    delete _temp_f;
    delete _ref_temp_f;
    delete _alpha_f;
    delete _velocity_f;
    delete _mach_f;
    delete _rho_air_f;
    delete _gamma_air_f;
    delete _dwdx_f;
    delete _dwdt_f;
    
    delete _thy;
    delete _thz;
    delete _hy_off;
    delete _hz_off;
    delete _E;
    delete _nu;
    delete _zero;
    delete _temp;
    delete _alpha;
    delete _velocity;
    delete _mach;
    delete _rho_air;
    delete _gamma_air;
    delete _dwdx;
    delete _dwdt;
    

    delete _eq_sys;
    delete _mesh;
    
    delete _discipline;
    delete _structural_sys;
    
    
}


