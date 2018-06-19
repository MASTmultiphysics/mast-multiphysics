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
#include "examples/fluid/cylinder_analysis_2D/cylinder_analysis_2d.h"
#include "examples/fluid/meshing/cylinder.h"
#include "examples/base/input_wrapper.h"


// libMesh includes
#include "libmesh/parallel_mesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/enum_elem_type.h"


MAST::Examples::CylinderAnalysis2D::
CylinderAnalysis2D(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::FluidExampleBase(comm_in) {
    
}


MAST::Examples::CylinderAnalysis2D::~CylinderAnalysis2D() {
    
}


void
MAST::Examples::CylinderAnalysis2D::_init_mesh() {

    _mesh              = new libMesh::ParallelMesh(this->comm());
    const unsigned int
    radial_divs         = (*_input)(_prefix+"n_radial_elems", "number of elements in the radial direction from cylinder to far-field", 20),
    quarter_divs        = (*_input)(_prefix+"n_quarter_elems", "number of elements in the quarter arc along the circumferencial direction", 20);

    const Real
    r                   = (*_input)("radius", "radius of the cylinder", 0.1),
    l_by_r              = (*_input)("l_by_r", "far-field distance to cylinder radius ratio", 5.),
    h_ff_by_h_r         = (*_input)("h_far_field_by_h_r", "relative element size at far-field boundary to element size at cylinder", 50.);
    
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

    // initialize the mesh
    MAST::Examples::CylinderMesh2D().mesh(r,
                                          r*l_by_r,
                                          radial_divs,
                                          quarter_divs,
                                          h_ff_by_h_r,
                                          *_mesh,
                                          e_type);
}



void
MAST::Examples::CylinderAnalysis2D::_init_loads() {
    
    bool
    if_viscous =
    (*_input)(_prefix+"if_viscous", "if the flow analysis should include viscosity", false);

    if (!if_viscous) {
    
        std::vector<unsigned int>
        slip      =  {3},
        no_slip,
        symm,
        far_field =  {1};
        _init_boundary_conditions(slip, no_slip, symm, far_field);
    }
    else {
        
        std::vector<unsigned int>
        slip,
        no_slip = {3},
        symm,
        far_field =  {1};
        _init_boundary_conditions(slip, no_slip, symm, far_field);
    }
}


