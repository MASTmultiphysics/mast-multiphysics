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
#include "elasticity/level_set_stress_assembly.h"
#include "elasticity/stress_output_base.h"
#include "elasticity/structural_system_initialization.h"
#include "base/nonlinear_system.h"
#include "level_set/level_set_intersection.h"


// libMesh includes
#include "libmesh/system.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"

MAST::LevelSetStressAssembly::LevelSetStressAssembly():
MAST::StressAssembly(),
_level_set     (nullptr),
_intersection  (nullptr) {
    
}


MAST::LevelSetStressAssembly::~LevelSetStressAssembly() {
    
}




void
MAST::LevelSetStressAssembly::
set_level_set_function(MAST::FieldFunction<Real>& level_set) {
    
    libmesh_assert(!_level_set);
    libmesh_assert(!_intersection);
    libmesh_assert(_system);
    
    _level_set    = &level_set;
    _intersection = new MAST::LevelSetIntersection(_system->system().get_mesh().max_elem_id(),
                                                   _system->system().get_mesh().max_node_id());
}



void
MAST::LevelSetStressAssembly::clear_level_set_function() {
    
    _level_set = nullptr;
    
    if (_intersection) {
        delete _intersection;
        _intersection = nullptr;
    }
}



extern void
get_max_stress_strain_values(const std::vector<MAST::StressStrainOutputBase::Data*>& data,
                             RealVectorX&           max_strain,
                             RealVectorX&           max_stress,
                             Real&                  max_vm,
                             const MAST::FunctionBase* p);




void MAST::LevelSetStressAssembly::
update_stress_strain_data(MAST::StressStrainOutputBase&       ops,
                          const libMesh::NumericVector<Real>& X) {
    
    // make sure it has been initialized
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(!_elem_ops);
    
    this->set_elem_operation_object(ops);
    
    RealVectorX sol;
    
    const MAST::NonlinearSystem&
    structural_sys = _system->system();
    
    libMesh::System&
    stress_sys = dynamic_cast<MAST::StructuralSystemInitialization*>(_system)->get_stress_sys();
    const std::vector<unsigned int>&
    stress_vars = dynamic_cast<MAST::StructuralSystemInitialization*>(_system)->get_stress_var_ids();
    
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = structural_sys.get_dof_map();
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(this->build_localized_vector(structural_sys,
                                                          X).release());
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    structural_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    structural_sys.get_mesh().active_local_elements_end();
    
    libMesh::dof_id_type
    dof_id   = 0;
    
    unsigned int
    sys_num =  stress_sys.number();
    
    RealVectorX
    max_strain_vals  =   RealVectorX::Zero(6),
    max_stress_vals  =   RealVectorX::Zero(6);
    
    Real
    max_vm_stress    =   0.;
    
    for ( ; el != end_el; el++) {
        
        const libMesh::Elem* elem = *el;
        
        _intersection->init(*_level_set, *elem, structural_sys.time);
        
        // process only if the element has any positive region
        if (!_intersection->if_elem_on_positive_phi()) continue;
        
        dof_map.dof_indices (elem, dof_indices);
        
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        // clear before calculating the data
        ops.clear();
        ops.init(*elem);
        ops.set_elem_solution(sol);
        ops.evaluate();
        ops.clear_elem();
        
        // get the stress-strain data map from the object
        const std::map<const libMesh::dof_id_type,
        std::vector<MAST::StressStrainOutputBase::Data*> >& output_map =
        ops.get_stress_strain_data();
        
        // make sure that only one element has been added to this data,
        // and that the element id is the same as the one being computed
        libmesh_assert_equal_to(output_map.size(), 1);
        libmesh_assert_equal_to(output_map.begin()->first, elem->id());
        
        // now iterate over all the elements and set the value in the
        // new system used for output
        std::map<const libMesh::dof_id_type,
        std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
        e_it    =  output_map.begin(),
        e_end   =  output_map.end();
        
        for ( ; e_it != e_end; e_it++) {
            
            get_max_stress_strain_values(e_it->second,
                                         max_strain_vals,
                                         max_stress_vals,
                                         max_vm_stress,
                                         nullptr);
            
            // set the values in the system
            // stress value
            dof_id     =   elem->dof_number(sys_num, stress_vars[12], 0);
            stress_sys.solution->set(dof_id, max_vm_stress);
            
            for (unsigned int i=0; i<6; i++) {
                // strain value
                dof_id     =   elem->dof_number(sys_num, stress_vars[i], 0);
                stress_sys.solution->set(dof_id, max_strain_vals(i));
                
                // stress value
                dof_id     =   elem->dof_number(sys_num, stress_vars[i+6], 0);
                stress_sys.solution->set(dof_id, max_stress_vals(i));
            }
        }
    }
    
    stress_sys.solution->close();
    this->clear_elem_operation_object();
}


