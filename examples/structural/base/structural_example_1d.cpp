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
#include "examples/structural/base/structural_example_1d.h"
#include "examples/base/multilinear_interpolation.h"
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

extern libMesh::LibMeshInit* __init;


MAST::Examples::StructuralExample1D::StructuralExample1D():
MAST::Examples::StructuralExampleBase() {

}


MAST::Examples::StructuralExample1D::~StructuralExample1D() {
 
}
        


void
MAST::Examples::StructuralExample1D::_init_mesh() {
    
    _mesh = new libMesh::SerialMesh(__init->comm());
    
    // identify the element type from the input file or from the order
    // of the element
    
    unsigned int
    nx_divs = (*_input)("nx_divs", 10);

    Real
    length  = (*_input)("length", 0.3);
    
    std::string
    t = (*_input)("elem_type", "edge2");

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
MAST::Examples::StructuralExample1D::_init_dirichlet_conditions() {

    this->_init_boundary_dirichlet_constraint(0, "left_constraint");
    this->_init_boundary_dirichlet_constraint(1, "right_constraint");
    
    _discipline->init_system_dirichlet_bc(*_sys);
}




void
MAST::Examples::StructuralExample1D::_init_section_property() {
    
    // the default behavior is without section offset
    this->_init_section_property_without_offset();
}

    
void
MAST::Examples::StructuralExample1D::_init_section_property_without_offset() {

    Real
    thy_v     =  (*_input)("thy",   0.06),
    thz_v     =  (*_input)("thz",   0.02);
    
    MAST::Parameter
    *thy      = new MAST::Parameter("thy", thy_v),
    *thz      = new MAST::Parameter("thz", thz_v),
    &zero     = this->get_parameter("zero");
    
    MAST::ConstantFieldFunction
    *thy_f    = new MAST::ConstantFieldFunction("hy",      *thy),
    *thz_f    = new MAST::ConstantFieldFunction("hz",      *thz),
    *hyoff_f  = new MAST::ConstantFieldFunction("hy_off",  zero),
    *hzoff_f  = new MAST::ConstantFieldFunction("hz_off",  zero);
    
    this->add_parameter(*thy);
    this->add_parameter(*thz);
    this->register_field_function(*thy_f);
    this->register_field_function(*thz_f);
    this->register_field_function(*hyoff_f);
    this->register_field_function(*hzoff_f);
    
    MAST::Solid1DSectionElementPropertyCard
    *p_card = new MAST::Solid1DSectionElementPropertyCard;
    
    libMesh::Point orientation;
    orientation(1) = 1.;
    p_card->y_vector() = orientation;
    _p_card = p_card;

    // set nonlinear strain if requested
    bool
    nonlinear  = (*_input)("if_nonlinear", false);
    if (nonlinear) p_card->set_strain(MAST::VON_KARMAN_STRAIN);
    
    p_card->add(*thy_f);
    p_card->add(*thz_f);
    p_card->add(*hyoff_f);
    p_card->add(*hzoff_f);
    p_card->set_material(*_m_card);
    p_card->init();
    _discipline->set_property_for_subdomain(0, *p_card);
}



void
MAST::Examples::StructuralExample1D::_init_section_property_with_offset() {
    
    Real
    thy_v     =  (*_input)("thy",   0.06),
    thz_v     =  (*_input)("thz",   0.02);
    
    MAST::Parameter
    *thy      = new MAST::Parameter("thy", thy_v),
    *thz      = new MAST::Parameter("thz", thz_v),
    &zero     = this->get_parameter("zero");
    
    MAST::ConstantFieldFunction
    *thy_f    = new MAST::ConstantFieldFunction("hy",      *thy),
    *thz_f    = new MAST::ConstantFieldFunction("hz",      *thz),
    *hzoff_f  = new MAST::ConstantFieldFunction("hz_off",  zero);
    
    MAST::SectionOffset
    *hyoff_f         = new MAST::SectionOffset("hy_off", *thy_f, 1.);
    
    this->add_parameter(*thy);
    this->add_parameter(*thz);
    this->register_field_function(*thy_f);
    this->register_field_function(*thz_f);
    this->register_field_function(*hyoff_f);
    this->register_field_function(*hzoff_f);
    
    MAST::Solid1DSectionElementPropertyCard
    *p_card = new MAST::Solid1DSectionElementPropertyCard;
    
    libMesh::Point orientation;
    orientation(1) = 1.;
    p_card->y_vector() = orientation;
    _p_card = p_card;
    
    // set nonlinear strain if requested
    bool
    nonlinear  = (*_input)("if_nonlinear", false);
    if (nonlinear) p_card->set_strain(MAST::VON_KARMAN_STRAIN);

    p_card->add(*thy_f);
    p_card->add(*thz_f);
    p_card->add(*hyoff_f);
    p_card->add(*hzoff_f);
    p_card->set_material(*_m_card);
    p_card->init();
    _discipline->set_property_for_subdomain(0, *p_card);
}

