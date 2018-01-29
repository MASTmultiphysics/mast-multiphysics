/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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
#include "base/parameter.h"
#include "base/boundary_condition_base.h"

// libMesh includes
#include "libmesh/string_to_enum.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"


MAST::Examples::ExampleBase::ExampleBase():
_initialized    (false),
_input          (nullptr) {
    
    MAST::Parameter
    *p = new MAST::Parameter("zero",  0.);
    this->add_parameter(*p);
}


MAST::Examples::ExampleBase::~ExampleBase() {
    
    {
        std::set<MAST::BoundaryConditionBase*>::iterator
        it   = _loadings.begin(),
        end  = _loadings.end();
        for ( ; it!=end; it++)
            delete *it;
    }

    {
        std::map<std::string, MAST::FunctionBase*>::iterator
        it   = _field_functions.begin(),
        end  = _field_functions.end();
        for ( ; it!=end; it++)
            delete it->second;
    }
    
    {
        std::map<std::string, MAST::Parameter*>::iterator
        it   = _parameters.begin(),
        end  = _parameters.end();
        for ( ; it!=end; it++)
            delete it->second;
    }
}



void
MAST::Examples::ExampleBase::init(GetPot& input) {
    
    libmesh_assert(!_initialized);
    
    _input       = &input;
    
    // FEType to initialize the system
    // get the order and type of element
    std::string
    order_str   = input( "fe_order",    "first"),
    family_str  = input("fe_family", "lagrange");
    
    libMesh::Order
    o  = libMesh::Utility::string_to_enum<libMesh::Order>(order_str);
    libMesh::FEFamily
    fe = libMesh::Utility::string_to_enum<libMesh::FEFamily>(family_str);
    _fetype = libMesh::FEType(o, fe);
    
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
    
    
    
    std::map<std::string, MAST::FunctionBase*>::iterator
    it   =  _field_functions.find(p.name()),
    end  =  _field_functions.end();
    
    // if the param was not found, then print the message
    if (it != end)
        libMesh::out
        << std::endl
        << "Function already exists by name: " << p.name() << std::endl;
    
    _field_functions.insert(std::map<std::string, MAST::FunctionBase*>::value_type(p.name(), &p));
}


void
MAST::Examples::ExampleBase::register_loading(MAST::BoundaryConditionBase& p) {
    
    
    // make sure that this already does not exist
    if (_loadings.count(&p))
        libMesh::out
        << std::endl
        << "Loading already exists." << std::endl;
    
    _loadings.insert(&p);
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


