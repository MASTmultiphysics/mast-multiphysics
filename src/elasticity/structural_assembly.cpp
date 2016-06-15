/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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
#include "elasticity/structural_assembly.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "base/physics_discipline_base.h"
#include "base/system_initialization.h"
#include "boundary_condition/point_load_condition.h"


// libMesh includes
#include "libmesh/petsc_vector.h"
#include "libmesh/implicit_system.h"
#include "libmesh/dof_map.h"


PetscErrorCode
MAST::_snes_structural_nonlinear_assembly_monitor_function(SNES snes,
                                                           PetscInt its,
                                                           PetscReal norm2,
                                                           void* ctx) {
    
    // convert the context pointer to the assembly object pointer
    MAST::StructuralNonlinearAssembly* assembly =
    (MAST::StructuralNonlinearAssembly*)ctx;
    
    // system for which this is being processed
    libMesh::System& sys = assembly->system();
    
    
    // if this is the first iteration, create a vector in system
    // as the initial approximation of the system solution
    if (its == 0) {
        sys.add_vector("old_solution");
    }
    else  {
        // use the previous solution as the base solution
        
        // get the solution update from SNES
        Vec dx;
        PetscErrorCode ierr = SNESGetSolutionUpdate(snes, &dx);
        libmesh_assert(!ierr);
        
        // attach the vector to a NumericVector object
        std::auto_ptr<libMesh::NumericVector<Real> >
        vec(new libMesh::PetscVector<Real>(dx, sys.comm())),
        vec_scaled(vec);
        
        // the solution update provided by SNES is -ve of the dX
        // hence, scale the vector by -1
        vec_scaled->scale(-1.);
        vec_scaled->close();
        
        // ask the assembly to update the incompatible mode solution
        // for all the 3D elements based on the solution update here
        assembly->update_incompatible_solution(sys.get_vector("old_solution"),
                                               *vec_scaled);
    }
    
    
    
    // finally, update the solution
    sys.get_vector("old_solution") = *sys.solution;
    sys.get_vector("old_solution").close();
    
    
    return 0;
}




MAST::StructuralAssembly::StructuralAssembly() {
    
}


MAST::StructuralAssembly::~StructuralAssembly() {
    
}


void
MAST::StructuralAssembly::
_assemble_point_loads(MAST::PhysicsDisciplineBase& discipline,
                      MAST::SystemInitialization& system,
                      libMesh::NumericVector<Real>& res) {
    
    // get a reference to the system, mesh and dof map
    libMesh::System& sys     = system.system();
    libMesh::DofMap& dof_map = sys.get_dof_map();
    
    // get a reference to the set of loads
    const MAST::PointLoadSetType& point_loads = discipline.point_loads();
    
    // iterate over the loads and process them
    MAST::PointLoadSetType::const_iterator
    load_it   = point_loads.begin(),
    load_end  = point_loads.end();
    
    for ( ; load_it != load_end; load_it++) {
        const MAST::PointLoadCondition& load = **load_it;
        const std::set<libMesh::Node*>& nodes = load.get_nodes();
        
        std::set<libMesh::Node*>::const_iterator
        n_it  = nodes.begin(),
        n_end = nodes.end();
        
        for ( ; n_it != n_end; n_it++) {
            
            // iterate over the nodes and the variables on them
            
            // get the dof id for the var on node
            
            // if the dof is not constrained, then add the
            
            // load to the residual vector
        }
        
    }
    
    res.close();
}



