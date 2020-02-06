/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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
#include "examples/base/example_base.h"
#include "examples/base/input_wrapper.h"
#include "base/parameter.h"
#include "base/boundary_condition_base.h"
#include "base/physics_discipline_base.h"
#include "base/nonlinear_system.h"
#include "base/system_initialization.h"
#include "property_cards/element_property_card_base.h"
#include "property_cards/material_property_card_base.h"


// libMesh includes
#include "libmesh/string_to_enum.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"


MAST::Examples::ExampleBase::ExampleBase(const libMesh::Parallel::Communicator& comm_in):
libMesh::ParallelObject (comm_in),
_initialized            (false),
_prefix                 (),
_input                  (nullptr),
_mesh                   (nullptr),
_eq_sys                 (nullptr),
_sys                    (nullptr),
_sys_init               (nullptr),
_discipline             (nullptr),
_m_card                 (nullptr),
_p_card                 (nullptr) {
    
    MAST::Parameter
    *p = new MAST::Parameter("zero",  0.);
    this->add_parameter(*p);
}


MAST::Examples::ExampleBase::~ExampleBase() {
    
    {
        std::set<MAST::BoundaryConditionBase*>::iterator
        it   = _boundary_conditions.begin(),
        end  = _boundary_conditions.end();
        for ( ; it!=end; it++)
            delete *it;
    }

    {
        std::set<MAST::FunctionBase*>::iterator
        it   = _field_functions.begin(),
        end  = _field_functions.end();
        for ( ; it!=end; it++)
            delete *it;
    }
    
    {
        std::map<std::string, MAST::Parameter*>::iterator
        it   = _parameters.begin(),
        end  = _parameters.end();
        for ( ; it!=end; it++)
            delete it->second;
    }
    
    if (!_initialized)
        return;
    
    delete _m_card;
    delete _p_card;
    
    delete _eq_sys;
    delete _mesh;
    
    delete _discipline;
    delete _sys_init;
}



void
MAST::Examples::ExampleBase::init(MAST::Examples::GetPotWrapper& input,
                                  const std::string& prefix) {
    
    libmesh_assert(!_initialized);
    
    _prefix      = prefix;
    _input       = &input;
    
    // call the initialization routines for each component
    _init_fetype();
    _init_mesh();
    _init_system_and_discipline();
    _init_dirichlet_conditions();
    _init_eq_sys();
    _init_material();
    _init_loads();
    _init_section_property();

    _initialized = true;
}



void
MAST::Examples::ExampleBase::add_parameter(MAST::Parameter& p) {
    
    
    
    std::map<std::string, MAST::Parameter*>::iterator
    it   =  _parameters.find(p.name()),
    end  =  _parameters.end();

    // if the param was not found, then print the message
    if (it != end)
        libMesh::out
        << std::endl
        << "Parameter already exists by name: " << p.name() << std::endl;
    
    _parameters.insert(std::map<std::string, MAST::Parameter*>::value_type(p.name(), &p));
}



MAST::Parameter&
MAST::Examples::ExampleBase::get_parameter(const std::string &nm) {
    
    std::map<std::string, MAST::Parameter*>::iterator
    it   =  _parameters.find(nm),
    end  =  _parameters.end();
    
    // if the param was not found, then print the message
    if (it == end)
        libMesh::out
        << std::endl
        << "Parameter not found by name: " << nm << std::endl;

    return *it->second;
}


void
MAST::Examples::ExampleBase::register_field_function(MAST::FunctionBase& p) {
    
    
    
    std::set<MAST::FunctionBase*>::iterator
    it   =  _field_functions.find(&p),
    end  =  _field_functions.end();
    
    // if the param was not found, then print the message
    if (it != end)
        libMesh::out
        << std::endl
        << "Function already exists: " << p.name() << std::endl;
    
    _field_functions.insert(&p);
}


MAST::FunctionBase&
MAST::Examples::ExampleBase::get_field_function(const std::string& nm) {

    std::set<MAST::FunctionBase*>::iterator
    it   =  _field_functions.begin(),
    end  =  _field_functions.end();

    for ( ; it != end; it++) {
        if ((*it)->name() == nm)
            break;
    }

    if (it == end) {
        // if it gets here, then we did not find the function
        libMesh::out
        << std::endl
        << "Function does not exits: " << nm << std::endl;
        libmesh_error();
    }
    
    return **it;
}



bool
MAST::Examples::ExampleBase::has_field_function(const std::string& nm) const {
    
    std::set<MAST::FunctionBase*>::iterator
    it   =  _field_functions.begin(),
    end  =  _field_functions.end();
    
    for ( ; it != end; it++) {
        if ((*it)->name() == nm)
            break;
    }
    
    if (it == end)
        return false;
    else
        return true;
}



void
MAST::Examples::ExampleBase::register_loading(MAST::BoundaryConditionBase& p) {
    
    
    // make sure that this already does not exist
    if (_boundary_conditions.count(&p))
        libMesh::out
        << std::endl
        << "Loading already exists." << std::endl;
    
    _boundary_conditions.insert(&p);
}


void
MAST::Examples::ExampleBase::add_load_parameter(MAST::Parameter& p) {
    
    std::map<MAST::Parameter*, const Real>::iterator
    it  = _load_parameters.find(&p),
    end = _load_parameters.end();
    
    if (it != end)
        libMesh::out
        << std::endl
        << "Parameter already exists by name: " << p.name() << std::endl;
    
    _load_parameters.insert(std::map<MAST::Parameter*, const Real>::value_type(&p, p()));
}



void
MAST::Examples::ExampleBase::update_load_parameters(Real scale) {
    
    libmesh_assert_greater_equal(scale, 0.);
    libmesh_assert_less_equal   (scale, 1.);
    
    std::map<MAST::Parameter*, const Real>::iterator
    it  = _load_parameters.begin(),
    end = _load_parameters.end();
    
    for ( ; it != end; it++)
        (*it->first) = scale*it->second;
}



void
MAST::Examples::ExampleBase::register_paramter_for_sensitivity(MAST::Parameter& p) {
    
    // make sure that the parameter does not already exist
    for (unsigned int i=0; i<_params_for_sensitivity.size(); i++)
        libmesh_assert_not_equal_to(_params_for_sensitivity[i], &p);
    
    _params_for_sensitivity.push_back(&p);
}


void
MAST::Examples::ExampleBase::_init_fetype() {
    
    // FEType to initialize the system
    // get the order and type of element
    std::string
    order_str   = (*_input)(_prefix+ "fe_order", "order of finite element shape basis functions",    "first"),
    family_str  = (*_input)(_prefix+"fe_family",      "family of finite element shape functions", "lagrange");
    
    libMesh::Order
    o  = libMesh::Utility::string_to_enum<libMesh::Order>(order_str);
    libMesh::FEFamily
    fe = libMesh::Utility::string_to_enum<libMesh::FEFamily>(family_str);
    _fetype = libMesh::FEType(o, fe);
}


