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
#include "level_set/level_set_void_solution.h"
#include "level_set/interface_dof_handler.h"
#include "level_set/level_set_intersection.h"
#include "base/physics_discipline_base.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "base/nonlinear_implicit_assembly_elem_operations.h"
#include "base/elem_base.h"

// libMesh includes
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/dof_map.h"


// monitor function for PETSc solver so that
// the incompatible solution can be updated after each converged iterate
PetscErrorCode
_snes_level_set_void_solution_assembly_monitor_function(SNES snes,
                                                        PetscInt its,
                                                        PetscReal norm2,
                                                        void* ctx) {
    
    // convert the context pointer to the assembly object pointer
    MAST::LevelSetVoidSolution* sol_update = (MAST::LevelSetVoidSolution*)ctx;
    
    // system for which this is being processed
    MAST::NonlinearSystem& sys = sol_update->get_assembly().system();
    
    if (sys.operation() != MAST::NonlinearSystem::NONLINEAR_SOLVE)
        return 0;
    
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
        sol_update->update_void_solution(sys.get_vector("old_solution"),
                                         *vec_scaled);
    }
    
    
    
    // finally, update the solution
    sys.get_vector("old_solution") = *sys.solution;
    sys.get_vector("old_solution").close();
    
    
    return 0;
}



MAST::LevelSetVoidSolution::LevelSetVoidSolution():
MAST::AssemblyBase::SolverMonitor(),
_assembly      (nullptr),
_intersection  (nullptr),
_dof_handler   (nullptr) {
    
    
}



MAST::LevelSetVoidSolution::~LevelSetVoidSolution() {
    
    this->clear();
}




void
MAST::LevelSetVoidSolution::init(MAST::AssemblyBase& assembly,
                                 MAST::LevelSetIntersection &intersection,
                                 MAST::LevelSetInterfaceDofHandler  &dof_handler) {

    // should be cleared before initialization
    libmesh_assert(!_assembly);
    
    _assembly    = &assembly;
    _intersection= &intersection;
    _dof_handler = &dof_handler;
    
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
                   _snes_level_set_void_solution_assembly_monitor_function,
                   (void*)this,
                   PETSC_NULL);
    
    libmesh_assert(!ierr);
}



void
MAST::LevelSetVoidSolution::clear() {
    
    // next, remove the monitor function from the snes object
    libMesh::PetscNonlinearSolver<Real> &petsc_nonlinear_solver =
    *(dynamic_cast<libMesh::PetscNonlinearSolver<Real>*>
      (dynamic_cast<MAST::NonlinearSystem&>(_assembly->system()).nonlinear_solver.get()));
    
    // get the SNES object
    SNES snes = petsc_nonlinear_solver.snes();
    
    PetscErrorCode ierr =
    SNESMonitorCancel(snes);
    libmesh_assert(!ierr);
    
    _assembly     = nullptr;
    _intersection = nullptr;
    _dof_handler  = nullptr;
}



void
MAST::LevelSetVoidSolution::
update_void_solution(libMesh::NumericVector<Real>& X,
                     libMesh::NumericVector<Real>& dX) {
    
    MAST::NonlinearImplicitAssemblyElemOperations
    &elem_ops  = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(_assembly->get_elem_ops());
    
    MAST::NonlinearSystem&
    sys        = _assembly->system();
    
    MAST::FieldFunction<Real>&
    phi        = _dof_handler->get_level_set_function();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX
    sol,
    dsol,
    sub_elem_vec,
    res;
    
    RealMatrixX
    sub_elem_mat,
    jac;
    
    std::vector<libMesh::dof_id_type>
    dof_indices,
    system_dofs,
    free_dofs;
    
    const libMesh::DofMap& dof_map = sys.get_dof_map();
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution (_assembly->build_localized_vector (sys,  X).release()),
    localized_dsolution(_assembly->build_localized_vector (sys, dX).release());


    libMesh::MeshBase::const_element_iterator       el     =
    sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    sys.get_mesh().active_local_elements_end();

    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem         = *el;
        if (_dof_handler->if_factor_element(*elem)) {
            
            dof_map.dof_indices (elem, dof_indices);
            
            // get the solution from the system vector
            unsigned int ndofs = (unsigned int)dof_indices.size();
            sol.setZero(ndofs);
            dsol.setZero(ndofs);
            sub_elem_vec.setZero(ndofs);
            res.setZero(ndofs);
            sub_elem_mat.setZero(ndofs, ndofs);
            jac.setZero(ndofs, ndofs);
            
            for (unsigned int i=0; i<dof_indices.size(); i++) {
                sol(i)   = (*localized_solution)(dof_indices[i]);
                dsol(i)  = (*localized_dsolution)(dof_indices[i]);
            }
            
            _dof_handler->solution_of_factored_element(*elem, sol);

            // get the intersection and compute the residual and jacobian
            // with contribution from all elements.
            _intersection->init(phi, *elem, sys.time);

            // the Jacobian is based on the homogenization method
            elem_ops.init(*elem);
            elem_ops.set_elem_solution(sol);
            elem_ops.elem_calculations(true, res, jac);
            elem_ops.clear_elem();
            jac *= _intersection->get_positive_phi_volume_fraction();
            res *= _intersection->get_positive_phi_volume_fraction();
            //res.setZero();

            
            /*const std::vector<const libMesh::Elem *> &
            elems_hi = _intersection->get_sub_elems_positive_phi();
            
            std::vector<const libMesh::Elem*>::const_iterator
            hi_sub_elem_it  = elems_hi.begin(),
            hi_sub_elem_end = elems_hi.end();
            
            for (; hi_sub_elem_it != hi_sub_elem_end; hi_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *hi_sub_elem_it;
                
                elem_ops.init(*sub_elem);
                elem_ops.set_elem_solution(sol);
                
                // if the element has been marked for factorization,
                // get the factorized jacobian and residual contributions
                elem_ops.elem_calculations(true, sub_elem_vec, sub_elem_mat);
                
                //jac += sub_elem_mat;
                res += sub_elem_vec;
                
                elem_ops.clear_elem();
            }*/

            _intersection->clear();
            
            _dof_handler->update_factored_element_solution(*elem,
                                                           res,
                                                           jac,
                                                           sol,
                                                           dsol,
                                                           sol);
        }
    }
}

