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
#include "examples/structural/beam_piston_theory_time_accurate/beam_piston_theory_time_accurate.h"
#include "base/nonlinear_system.h"

// libMesh includes
#include "libmesh/unstructured_mesh.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"


MAST::Examples::BeamPistonTheoryTimeAccurateAnalysis::BeamPistonTheoryTimeAccurateAnalysis():
MAST::Examples::StructuralExample1D() {
    
}





MAST::Examples::BeamPistonTheoryTimeAccurateAnalysis::~BeamPistonTheoryTimeAccurateAnalysis() {
    
}



void
MAST::Examples::BeamPistonTheoryTimeAccurateAnalysis::_init_loads() {

    this->_init_temperature_load();
    this->_init_piston_theory_load();
}




void
MAST::Examples::BeamPistonTheoryTimeAccurateAnalysis::initialize_solution() {
    
    libmesh_assert(_initialized);
        
    // initialize the solution before solving
    libMesh::MeshBase::node_iterator
    n_id  = _mesh->nodes_begin(),
    n_end = _mesh->nodes_end();
    unsigned int
    dof = 0;

    Real
    x      = 0.,
    length = (*_input)("length", 0.3),
    pi     = acos(-1.);
    
    for ( ; n_id != n_end; n_id++) {
        
        libMesh::Node& n = **n_id;
        x = n(0);
        
        // set v
        dof = n.dof_number(_sys->number(), 1, 0);
        _sys->solution->set(dof, 0.1*sin(pi*x/length));
        // set theta_z
        dof = n.dof_number(_sys->number(), 5, 0);
        _sys->solution->set(dof, 0.1*pi/length*cos(pi*x/length));
        _sys->solution->close();
    }
}



