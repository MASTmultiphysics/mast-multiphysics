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
#include "examples/structural/nonlinear_circular_cantilever/circular_cantilever.h"
#include "examples/structural/base/thermal_stress_jacobian_scaling_function.h"
#include "examples/base/multilinear_interpolation.h"
#include "examples/base/input_wrapper.h"
#include "base/nonlinear_system.h"
#include "base/constant_field_function.h"
#include "base/parameter.h"
#include "base/physics_discipline_base.h"
#include "property_cards/isotropic_element_property_card_3D.h"
#include "boundary_condition/dirichlet_boundary_condition.h"

// libMesh includes
#include "libmesh/serial_mesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"



MAST::Examples::StructuralExample3D::
StructuralExample3D(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::StructuralExampleBase(comm_in) {
    
}


MAST::Examples::StructuralExample3D::~StructuralExample3D() {
    
}



void
MAST::Examples::StructuralExample3D::_init_mesh() {
    
    _mesh = new libMesh::SerialMesh(this->comm());
    
    // identify the element type from the input file or from the order
    // of the element
    
    unsigned int
    nx_divs = (*_input)(_prefix+"nx_divs", "number of elements along x-axis", 10),
    ny_divs = (*_input)(_prefix+"ny_divs", "number of elements along x-axis", 5),
    nz_divs = (*_input)(_prefix+"nz_divs", "number of elements along y-axis", 5);
    
    Real
    length  = (*_input)(_prefix+ "length", "length of domain along x-axis", 0.3),
    width   = (*_input)(_prefix+  "width", "length of domain along y-axis", 0.005),
    height  = (*_input)(_prefix+ "height", "length of domain along y-axis", 0.01);
    
    std::string
    t = (*_input)(_prefix+"elem_type", "type of geometric element in the mesh", "hex8");
    
    libMesh::ElemType
    e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
    
    // if high order FE is used, libMesh requires atleast a second order
    // geometric element.
    if (_fetype.order > 1 && e_type == libMesh::HEX8)
        e_type = libMesh::HEX27;
    else if (_fetype.order > 1 && e_type == libMesh::TET4)
        e_type = libMesh::TET10;
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_cube(*_mesh,
                                               nx_divs, ny_divs, nz_divs,
                                               0, length,
                                               0, width,
                                               0, height,
                                               e_type);
}



void
MAST::Examples::StructuralExample3D::_init_dirichlet_conditions() {
    
    // side id associations based on libMesh's indexing for hex8 and
    // cube mesh generation
    //this->_init_boundary_dirichlet_constraint(0, "back_constraint");
    //this->_init_boundary_dirichlet_constraint(1, "bottom_constraint");
    //this->_init_boundary_dirichlet_constraint(2, "right_constraint");
    //this->_init_boundary_dirichlet_constraint(3, "top_constraint");
    this->_init_boundary_dirichlet_constraint(4, "left_constraint");
    //this->_init_boundary_dirichlet_constraint(5, "front_constraint");

    _discipline->init_system_dirichlet_bc(*_sys);
}




void
MAST::Examples::StructuralExample3D::_init_section_property() {
    
    MAST::IsotropicElementPropertyCard3D
    *p_card   = new MAST::IsotropicElementPropertyCard3D;
    
    _p_card   = p_card;
    
    // set nonlinear strain if requested
    bool
    nonlinear = (*_input)(_prefix+"if_nonlinear", "flag to turn on/off nonlinear strain", false);
    if (nonlinear) p_card->set_strain(MAST::NONLINEAR_STRAIN);
    
    p_card->set_material(*_m_card);
    _discipline->set_property_for_subdomain(0, *p_card);
}


void
MAST::Examples::StructuralExample3D::_init_loads() {
    
    _init_pressure_load(true, 3);
    //_init_temperature_load();
}
