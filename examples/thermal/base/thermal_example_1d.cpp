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
#include "examples/thermal/base/thermal_example_1d.h"
#include "examples/base/input_wrapper.h"
#include "base/nonlinear_system.h"
#include "base/constant_field_function.h"
#include "base/parameter.h"
#include "base/physics_discipline_base.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "boundary_condition/dirichlet_boundary_condition.h"

// libMesh includes
#include "libmesh/serial_mesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"



MAST::Examples::ThermalExample1D::
ThermalExample1D(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::ThermalExampleBase(comm_in) {
    
}


MAST::Examples::ThermalExample1D::~ThermalExample1D() {
    
}



void
MAST::Examples::ThermalExample1D::_init_mesh() {
    
    _mesh = new libMesh::SerialMesh(this->comm());
    
    // identify the element type from the input file or from the order
    // of the element
    
    unsigned int
    nx_divs = (*_input)(_prefix+"nx_divs", "number of elements along x-axis", 10);
    
    Real
    length  = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3);
    
    std::string
    t = (*_input)(_prefix+"elem_type", "type of geometric element in the mesh", "edge2");
    
    libMesh::ElemType
    e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
    
    // if high order FE is used, libMesh requires atleast a second order
    // geometric element.
    if (_fetype.order > 1 && e_type == libMesh::EDGE2)
        e_type = libMesh::EDGE3;
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_line(*_mesh, nx_divs, 0, length, e_type);
}



void
MAST::Examples::ThermalExample1D::_init_dirichlet_conditions() {
    
    this->_init_boundary_dirichlet_constraint(0, "left_constraint");
    this->_init_boundary_dirichlet_constraint(1, "right_constraint");
    
    _discipline->init_system_dirichlet_bc(*_sys);
}




void
MAST::Examples::ThermalExample1D::_init_section_property() {
    
    Real
    thy_v     =  (*_input)(_prefix+"thy", "thickness of element along y-axis",  0.06),
    thz_v     =  (*_input)(_prefix+"thz", "thickness of element along z-axis",  0.02);
    
    MAST::Parameter
    *thy      = new MAST::Parameter("thy", thy_v),
    *thz      = new MAST::Parameter("thz", thz_v);
    
    MAST::ConstantFieldFunction
    *thy_f    = new MAST::ConstantFieldFunction("hy",      *thy),
    *thz_f    = new MAST::ConstantFieldFunction("hz",      *thz);
    
    this->add_parameter(*thy);
    this->add_parameter(*thz);
    this->register_field_function(*thy_f);
    this->register_field_function(*thz_f);
    
    MAST::Solid1DSectionElementPropertyCard
    *p_card = new MAST::Solid1DSectionElementPropertyCard;
    
    libMesh::Point orientation;
    orientation(1) = 1.;
    p_card->y_vector() = orientation;
    _p_card = p_card;
    
    p_card->add(*thy_f);
    p_card->add(*thz_f);
    p_card->set_material(*_m_card);
    p_card->init();
    _discipline->set_property_for_subdomain(0, *p_card);
}


