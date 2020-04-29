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
#include "level_set/heaviside_elem_homogenization_function.h"
#include "level_set/level_set_elem_base.h"
#include "level_set/filter_base.h"
#include "level_set/level_set_parameter.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "mesh/geom_elem.h"

// libMesh includes
#include "libmesh/elem.h"


MAST::HeavisideElemHomogenizedDensityFunction::
HeavisideElemHomogenizedDensityFunction(const std::string& nm):
MAST::HomogenizedDensityFunctionBase(nm),
_width          (0.1) {
    
}



MAST::HeavisideElemHomogenizedDensityFunction::
~HeavisideElemHomogenizedDensityFunction() {
    
}


void
MAST::HeavisideElemHomogenizedDensityFunction::
initialize_element_volume_fractions() {
    
    libmesh_assert(_analysis_mesh);
    libmesh_assert(_elem_volume_fraction.empty());

    // make sure that the system uses Lagrange shape functions for
    // the level set.
    libmesh_assert(_level_set_sys->system().variable_type(0).family == libMesh::LAGRANGE);
    
    // Next, compute the value of level-set at the nodes of the analysis mesh
    // and obtain the value of level-set from the mesh function.
    // These values will be provided to the elmeent for computation of
    // the homogenized volume fraction.

    RealVectorX
    phi;
    
    unsigned int
    n_nodes = 0;
    
    Real
    v = 0.;
    
    libMesh::MeshBase::const_element_iterator
    e_it   =  _analysis_mesh->active_local_elements_begin(),
    e_end  =  _analysis_mesh->active_local_elements_end();
    
    for ( ; e_it != e_end; e_it++) {

        const libMesh::Elem* e = *e_it;

        n_nodes = e->n_nodes();
        phi.setZero(n_nodes);
        
        // get the nodal values for this element
        for (unsigned int i=0; i<n_nodes; i++) {
            (*_level_set)(*e->node_ptr(i), 0., v);
            phi(i) = v;
        }
        
        // compute the integral of Heaviside function to approximate the
        // volume fraction
        MAST::GeomElem             geom_elem;
        geom_elem.init(*e, *_level_set_sys);
        MAST::LevelSetElementBase  level_set_elem(*_level_set_sys, geom_elem);
        level_set_elem.set_solution(phi);
        v = level_set_elem.homogenized_volume_fraction(_width);
        
        // store the value for this elmeent
        _elem_volume_fraction[e] = v;
    }
}



void
MAST::HeavisideElemHomogenizedDensityFunction::
initialize_element_volume_fraction_sensitivity(const MAST::FunctionBase& f) {
    
    libmesh_assert(_analysis_mesh);
    libmesh_assert(_elem_volume_fraction_sensitivity[&f].empty());
    libmesh_assert(f.is_topology_parameter());
    
    // make sure that the system uses Lagrange shape functions for
    // the level set.
    libmesh_assert(_level_set_sys->system().variable_type(0).family == libMesh::LAGRANGE);
    
    // Next, compute the value of level-set at the nodes of the analysis mesh
    // and obtain the value of level-set from the mesh function.
    // These values will be provided to the elmeent for computation of
    // the homogenized volume fraction.

    RealVectorX
    phi,
    phi_sens;
    
    unsigned int
    n_nodes = 0;
    
    Real
    v = 0.;
    
    const MAST::LevelSetParameter
    &p_ls = dynamic_cast<const MAST::LevelSetParameter&>(f);
    
    libMesh::MeshBase::const_element_iterator
    e_it   =  _analysis_mesh->active_local_elements_begin(),
    e_end  =  _analysis_mesh->active_local_elements_end();
    
    for ( ; e_it != e_end; e_it++) {

        const libMesh::Elem* e = *e_it;

        if (_filter->if_elem_in_domain_of_influence(*e, *p_ls.level_set_node())) {
            
            n_nodes = e->n_nodes();
            phi.setZero(n_nodes);
            phi_sens.setZero(n_nodes);
            
            // get the nodal values for this element
            for (unsigned int i=0; i<n_nodes; i++) {
                
                (*_level_set)(*e->node_ptr(i), 0., v);
                phi(i) = v;
                _level_set->derivative(f, *e->node_ptr(i), 0., v);
                phi_sens(i) = v;
            }
            
            // compute the integral of Heaviside function to approximate the
            // volume fraction
            MAST::GeomElem             geom_elem;
            geom_elem.init(*e, *_level_set_sys);
            MAST::LevelSetElementBase  level_set_elem(*_level_set_sys, geom_elem);
            level_set_elem.set_solution(phi);
            level_set_elem.set_solution(phi_sens, true);
            v = level_set_elem.homogenized_volume_fraction_sensitivity(_width);
            
            // store the value for this elmeent
            _elem_volume_fraction_sensitivity[&f][e] = v;
        }
        else
            _elem_volume_fraction_sensitivity[&f][e] = 0.;
    }
}


