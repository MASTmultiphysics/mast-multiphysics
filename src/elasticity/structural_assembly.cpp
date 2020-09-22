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
#include "elasticity/structural_assembly.h"
#include "base/physics_discipline_base.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "base/assembly_elem_operation.h"
#include "elasticity/structural_element_base.h"
#include "mesh/geom_elem.h"

// libMesh includes
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/dof_map.h"


// monitor function for PETSc solver so that
// the incompatible solution can be updated after each converged iterate
PetscErrorCode
_snes_structural_nonlinear_assembly_monitor_function(SNES snes,
                                                     PetscInt its,
                                                     PetscReal norm2,
                                                     void* ctx) {
    
    // convert the context pointer to the assembly object pointer
    MAST::AssemblyBase* assembly = (MAST::AssemblyBase*)ctx;
    
    // system for which this is being processed
    MAST::NonlinearSystem& sys = assembly->system();

    MAST::StructuralAssembly*
    str_assembly = dynamic_cast<MAST::StructuralAssembly*>(assembly->get_solver_monitor());
    
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
        std::unique_ptr<libMesh::NumericVector<Real> >
        vec(new libMesh::PetscVector<Real>(dx, sys.comm())),
        vec_scaled(vec->clone().release());
        
        // the solution update provided by SNES is -ve of the dX
        // hence, scale the vector by -1
        vec_scaled->scale(-1.);
        vec_scaled->close();
        
        // ask the assembly to update the incompatible mode solution
        // for all the 3D elements based on the solution update here
        str_assembly->update_incompatible_solution(sys.get_vector("old_solution"),
                                                   *vec_scaled);
    }
    
    
    
    // finally, update the solution
    sys.get_vector("old_solution") = *sys.solution;
    sys.get_vector("old_solution").close();
    
    
    return 0;
}




MAST::StructuralAssembly::StructuralAssembly():
MAST::AssemblyBase::SolverMonitor(),
_assembly (nullptr) {
    
}


MAST::StructuralAssembly::~StructuralAssembly() {
    
}


void
MAST::StructuralAssembly::set_elem_incompatible_sol(MAST::StructuralElementBase &elem) {
    
    const libMesh::Elem& e = elem.elem().get_reference_elem();

    if (elem.if_incompatible_modes()) {
        
        // init the solution if it is not currently set
        if (!_incompatible_sol.count(&e))
            _incompatible_sol[&e] = RealVectorX::Zero(elem.incompatible_mode_size());
        
        // use the solution currently available 
        elem.set_incompatible_mode_solution(_incompatible_sol[&e]);
    }
}


void
MAST::StructuralAssembly::init(MAST::AssemblyBase& assembly) {

    // should be cleared before initialization
    libmesh_assert(!_assembly);
    
    _assembly = &assembly;
    
    // get the nonlinear solver SNES object from System and
    // add a monitor to it so that it can be used to update the
    // incompatible mode solution after each update
    
    libMesh::PetscNonlinearSolver<Real> &petsc_nonlinear_solver =
    *(dynamic_cast<libMesh::PetscNonlinearSolver<Real>*>
      (dynamic_cast<MAST::NonlinearSystem&>(assembly.system()).nonlinear_solver.get()));
    
    
    // initialize the solver before getting the snes object
    if (libMesh::on_command_line("--solver_system_names"))
        petsc_nonlinear_solver.init((assembly.system().name()+"_").c_str());
    
    // get the SNES object
    SNES snes = petsc_nonlinear_solver.snes();
    
    PetscErrorCode ierr =
    SNESMonitorSet(snes,
                   _snes_structural_nonlinear_assembly_monitor_function,
                   (void*)this,
                   PETSC_NULL);
    
    libmesh_assert(!ierr);
}



void
MAST::StructuralAssembly::clear() {
    
    // next, remove the monitor function from the snes object
    libMesh::PetscNonlinearSolver<Real> &petsc_nonlinear_solver =
    *(dynamic_cast<libMesh::PetscNonlinearSolver<Real>*>
      (dynamic_cast<MAST::NonlinearSystem&>(_assembly->system()).nonlinear_solver.get()));
    
    // get the SNES object
    SNES snes = petsc_nonlinear_solver.snes();
    
    PetscErrorCode ierr =
    SNESMonitorCancel(snes);
    libmesh_assert(!ierr);
    
    _assembly = nullptr;
}



void
MAST::StructuralAssembly::update_incompatible_solution(libMesh::NumericVector<Real>& X,
                                                       libMesh::NumericVector<Real>& dX) {
    
    
    // iterate over each element and ask the 3D elements to update
    // their local solutions
    
    MAST::NonlinearSystem& sys = _assembly->system();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol, dsol;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = sys.get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution(_assembly->build_localized_vector(sys,
                                                         X).release()),
    localized_dsolution(_assembly->build_localized_vector(sys,
                                                          dX).release());
    
    
    // if a solution function is attached, initialize it
    //if (_sol_function)
    //    _sol_function->init( X, false);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem         = *el;
        MAST::AssemblyElemOperations& ops = _assembly->get_elem_ops();
        
        MAST::GeomElem geom_elem;
        ops.set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, _assembly->system_init());
        
        ops.init(geom_elem);

        MAST::StructuralElementBase& p_elem =
        dynamic_cast<MAST::StructuralElementBase&>(ops.get_physics_elem());
        
        if (p_elem.if_incompatible_modes()) {
            
            dof_map.dof_indices (elem, dof_indices);
            
            // get the solution
            unsigned int ndofs = (unsigned int)dof_indices.size();
            sol.setZero(ndofs);
            dsol.setZero(ndofs);
            
            for (unsigned int i=0; i<dof_indices.size(); i++) {
                sol(i)  = (*localized_solution) (dof_indices[i]);
                dsol(i) = (*localized_dsolution)(dof_indices[i]);
            }
            
            p_elem.set_solution(sol);
            p_elem.set_incompatible_mode_solution(_incompatible_sol[elem]);
            
            //if (_sol_function)
            //    p_elem.attach_active_solution_function(*_sol_function);
            
            // perform the element level calculations
            p_elem.update_incompatible_mode_solution(dsol);
            
            p_elem.detach_active_solution_function();
        }
    }
    
    
    // if a solution function is attached, clear it
    //if (_sol_function)
    //    _sol_function->clear();
}


//void
//MAST::StructuralAssembly::
//_assemble_point_loads(MAST::PhysicsDisciplineBase& discipline,
//                      MAST::SystemInitialization& system,
//                      libMesh::NumericVector<Real>& res) {
//    
//    // get a reference to the system, mesh and dof map
//    MAST::NonlinearSystem& sys     = system.system();
//    libMesh::DofMap& dof_map       = sys.get_dof_map();
//    
//    // get a reference to the set of loads
//    const MAST::PointLoadSetType& point_loads = discipline.point_loads();
//    
//    // iterate over the loads and process them
//    MAST::PointLoadSetType::const_iterator
//    load_it   = point_loads.begin(),
//    load_end  = point_loads.end();
//    
//    for ( ; load_it != load_end; load_it++) {
//        const MAST::PointLoadCondition& load = **load_it;
//        const std::set<libMesh::Node*>& nodes = load.get_nodes();
//        
//        std::set<libMesh::Node*>::const_iterator
//        n_it  = nodes.begin(),
//        n_end = nodes.end();
//        
//        for ( ; n_it != n_end; n_it++) {
//            
//            // iterate over the nodes and the variables on them
//            
//            // get the dof id for the var on node
//            
//            // if the dof is not constrained, then add the
//            
//            // load to the residual vector
//        }
//        
//    }
//    
//    res.close();
//}



