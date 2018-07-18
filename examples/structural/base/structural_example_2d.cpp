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
#include "examples/structural/base/structural_example_2d.h"
#include "examples/structural/base/thermal_stress_jacobian_scaling_function.h"
#include "examples/base/multilinear_interpolation.h"
#include "examples/base/input_wrapper.h"
#include "base/nonlinear_system.h"
#include "base/constant_field_function.h"
#include "base/parameter.h"
#include "base/physics_discipline_base.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "boundary_condition/dirichlet_boundary_condition.h"

// libMesh includes
#include "libmesh/serial_mesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"



MAST::Examples::StructuralExample2D::
StructuralExample2D(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::StructuralExampleBase(comm_in) {
    
}


MAST::Examples::StructuralExample2D::~StructuralExample2D() {
    
}



void
MAST::Examples::StructuralExample2D::_init_mesh() {
    
    _mesh = new libMesh::SerialMesh(this->comm());
    
    // identify the element type from the input file or from the order
    // of the element
    
    unsigned int
    nx_divs = (*_input)(_prefix+"nx_divs", "number of elements along x-axis", 10),
    ny_divs = (*_input)(_prefix+"ny_divs", "number of elements along y-axis", 10);
    
    Real
    length  = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
    width   = (*_input)(_prefix+ "width", "length of domain along y-axis", 0.3);
    
    std::string
    t = (*_input)(_prefix+"elem_type", "type of geometric element in the mesh", "quad4");
    
    libMesh::ElemType
    e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
    
    // if high order FE is used, libMesh requires atleast a second order
    // geometric element.
    if (_fetype.order > 1 && e_type == libMesh::QUAD4)
        e_type = libMesh::QUAD9;
    else if (_fetype.order > 1 && e_type == libMesh::TRI3)
        e_type = libMesh::TRI6;

    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_square(*_mesh,
                                                 nx_divs, ny_divs,
                                                 0, length,
                                                 0, width,
                                                 e_type);
}



void
MAST::Examples::StructuralExample2D::_init_dirichlet_conditions() {

    this->_init_boundary_dirichlet_constraint(0, "bottom_constraint");
    this->_init_boundary_dirichlet_constraint(1, "right_constraint");
    this->_init_boundary_dirichlet_constraint(2, "top_constraint");
    this->_init_boundary_dirichlet_constraint(3, "left_constraint");

    _discipline->init_system_dirichlet_bc(*_sys);
}




void
MAST::Examples::StructuralExample2D::_init_section_property() {
    
    // the default behavior is without section offset
    this->_init_section_property_without_offset();
}


void
MAST::Examples::StructuralExample2D::_init_section_property_without_offset() {
    
    Real
    th_v      =  (*_input)(_prefix+"th", "thickness of 2D element",  0.001);
    
    MAST::Parameter
    *th       = new MAST::Parameter("th", th_v),
    &zero     = this->get_parameter("zero");
    
    MAST::ConstantFieldFunction
    *th_f     = new MAST::ConstantFieldFunction("h",       *th),
    *hoff_f   = new MAST::ConstantFieldFunction("off",    zero);
    
    
    this->add_parameter(*th);
    this->register_field_function(*th_f);
    this->register_field_function(*hoff_f);
    
    MAST::Solid2DSectionElementPropertyCard
    *p_card   = new MAST::Solid2DSectionElementPropertyCard;
    
    _p_card   = p_card;

    // set nonlinear strain if requested
    bool
    nonlinear = (*_input)(_prefix+"if_nonlinear", "flag to turn on/off nonlinear strain", false);
    if (nonlinear) p_card->set_strain(MAST::NONLINEAR_STRAIN);
    
    p_card->add(*th_f);
    p_card->add(*hoff_f);
    p_card->set_material(*_m_card);
    _discipline->set_property_for_subdomain(0, *p_card);
}



void
MAST::Examples::StructuralExample2D::_init_section_property_with_offset() {
    
    Real
    th_v      =  (*_input)(_prefix+"th", "thickness of 2D element",  0.006);

    MAST::Parameter
    *th       = new MAST::Parameter("th", th_v);
    
    MAST::ConstantFieldFunction
    *th_f     = new MAST::ConstantFieldFunction("h",       *th);

    MAST::SectionOffset
    *hoff_f   = new MAST::SectionOffset("off", *th_f, 1.);

    this->add_parameter(*th);
    this->register_field_function(*th_f);
    this->register_field_function(*hoff_f);
    
    MAST::Solid2DSectionElementPropertyCard
    *p_card   = new MAST::Solid2DSectionElementPropertyCard;
    
    _p_card   = p_card;
    
    // set nonlinear strain if requested
    bool
    nonlinear = (*_input)(_prefix+"if_nonlinear", "flag to turn on/off nonlinear strain", false);
    if (nonlinear) p_card->set_strain(MAST::NONLINEAR_STRAIN);
    
    p_card->add(*th_f);
    p_card->add(*hoff_f);
    p_card->set_material(*_m_card);
    _discipline->set_property_for_subdomain(0, *p_card);
}

