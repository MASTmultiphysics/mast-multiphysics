/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_assembly.h"
#include "property_cards/element_property_card_base.h"
#include "base/physics_discipline_base.h"
#include "base/system_initialization.h"
#include "base/mesh_field_function.h"
#include "numerics/utility.h"
#include "base/real_output_function.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/parameter_vector.h"



MAST::StructuralNonlinearAssembly::
StructuralNonlinearAssembly():
MAST::NonlinearImplicitAssembly() {
    
}



MAST::StructuralNonlinearAssembly::
~StructuralNonlinearAssembly() {
    
}



void
MAST::StructuralNonlinearAssembly::
residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                       libMesh::NumericVector<Real>* R,
                       libMesh::SparseMatrix<Real>*  J,
                       libMesh::NonlinearImplicitSystem& S) {
    
    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &(nonlin_sys));
    
    if (R) R->zero();
    if (J) J->zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(_build_localized_vector(nonlin_sys,
                                                     X).release());
    
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        physics_elem.reset(_build_elem(*elem).release());

        MAST::StructuralElementBase& p_elem =
        dynamic_cast<MAST::StructuralElementBase&>(*physics_elem);

        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        physics_elem->set_solution    (sol);
        physics_elem->set_velocity    (vec); // set to zero vector for a quasi-steady analysis
        physics_elem->set_acceleration(vec); // set to zero vector for a quasi-steady analysis
        
        
        // set the incompatible mode solution if required by the
        // element
        if (p_elem.if_incompatible_modes()) {
            // check if the vector exists in the map
            if (!_incompatible_sol.count(elem))
                _incompatible_sol[elem] = RealVectorX::Zero(p_elem.incompatible_mode_size());
            p_elem.set_incompatible_mode_solution(_incompatible_sol[elem]);
        }
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);
        
        //_check_element_numerical_jacobian(*physics_elem, sol);
        
        // perform the element level calculations
        _elem_calculations(*physics_elem,
                           J!=nullptr?true:false,
                           vec, mat);
        
        physics_elem->detach_active_solution_function();
        
        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        DenseRealMatrix m;
        if (R)
            MAST::copy(v, vec);
        if (J)
            MAST::copy(m, mat);
        
        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        if (R && J)
            dof_map.constrain_element_matrix_and_vector(m, v, dof_indices);
        else if (R)
            dof_map.constrain_element_vector(v, dof_indices);
        else
            dof_map.constrain_element_matrix(m, dof_indices);
        
        // add to the global matrices
        if (R) R->add_vector(v, dof_indices);
        if (J) J->add_matrix(m, dof_indices);
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    if (R) R->close();
    if (J) J->close();
}





bool
MAST::StructuralNonlinearAssembly::
sensitivity_assemble (const libMesh::ParameterVector& parameters,
                      const unsigned int i,
                      libMesh::NumericVector<Real>& sensitivity_rhs) {
    
    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    sensitivity_rhs.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(_build_localized_vector(nonlin_sys,
                                                     *nonlin_sys.solution).release());
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( *nonlin_sys.solution);
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        physics_elem.reset(_build_elem(*elem).release());
        
        MAST::StructuralElementBase& p_elem =
        dynamic_cast<MAST::StructuralElementBase&>(*physics_elem);

        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        physics_elem->sensitivity_param = _discipline->get_parameter(&(parameters[i].get()));
        physics_elem->set_solution    (sol);
        physics_elem->set_velocity    (vec); // set to zero vector for a quasi-steady analysis
        physics_elem->set_acceleration(vec); // set to zero vector for a quasi-steady analysis
        
        // set the incompatible mode solution if required by the
        // element
        if (p_elem.if_incompatible_modes()) {
            // check if the vector exists in the map
            if (!_incompatible_sol.count(elem))
                _incompatible_sol[elem] = RealVectorX::Zero(p_elem.incompatible_mode_size());
            p_elem.set_incompatible_mode_solution(_incompatible_sol[elem]);
        }

        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);
        
        // perform the element level calculations
        _elem_sensitivity_calculations(*physics_elem, false, vec, mat);
        
        // the sensitivity method provides sensitivity of the residual.
        // Hence, this is multiplied with -1 to make it the RHS of the
        // sensitivity equations.
        vec *= -1.;
        
        physics_elem->detach_active_solution_function();
        
        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        MAST::copy(v, vec);
        
        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        dof_map.constrain_element_vector(v, dof_indices);
        
        // add to the global matrices
        sensitivity_rhs.add_vector(v, dof_indices);
    }
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->clear();
    
    sensitivity_rhs.close();
    
    return true;
}





void
MAST::StructuralNonlinearAssembly::StructuralNonlinearAssembly::
update_incompatible_solution(libMesh::NumericVector<Real>& X,
                             libMesh::NumericVector<Real>& dX) {
    
    
    // iterate over each element and ask the 3D elements to update
    // their local solutions
    
    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol, dsol;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> >
    localized_solution(_build_localized_vector(nonlin_sys,
                                               X).release()),
    localized_dsolution(_build_localized_vector(nonlin_sys,
                                                dX).release());
    
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        physics_elem.reset(_build_elem(*elem).release());
        
        MAST::StructuralElementBase& p_elem =
        dynamic_cast<MAST::StructuralElementBase&>(*physics_elem);

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
            
            if (_sol_function)
                p_elem.attach_active_solution_function(*_sol_function);
            
            // perform the element level calculations
            p_elem.update_incompatible_mode_solution(dsol);
            
            p_elem.detach_active_solution_function();
        }
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
}




void
MAST::StructuralNonlinearAssembly::
calculate_outputs(const libMesh::NumericVector<Real>& X) {

    
    // look through the outputs for structural compliance. If present,
    // then evaluate it separately
    MAST::VolumeOutputMapType::iterator
    it    = _discipline->volume_output().begin(),
    end   = _discipline->volume_output().end();
    
    
    for ( ; it != end; it++)
        if (it->second->type() == MAST::STRUCTURAL_COMPLIANCE) {
            
            libMesh::NonlinearImplicitSystem& sys =
            dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
            
            MAST::RealOutputFunction& output =
            dynamic_cast<MAST::RealOutputFunction&>(*it->second);
            
            _calculate_compliance(X, sys, output);
        }
    
    // now call the parent's method for calculation of the other outputs
    MAST::NonlinearImplicitAssembly::calculate_outputs(X);
    
}






void
MAST::StructuralNonlinearAssembly::
calculate_output_sensitivity(libMesh::ParameterVector& params,
                             const bool if_total_sensitivity,
                             const libMesh::NumericVector<Real>& X) {


    // look through the outputs for structural compliance. If present,
    // then evaluate it separately
    MAST::VolumeOutputMapType::iterator
    it    = _discipline->volume_output().begin(),
    end   = _discipline->volume_output().end();
    
    std::auto_ptr<libMesh::NumericVector<Real> >
    localized_solution_sensitivity;
    
    for ( ; it != end; it++)
        if (it->second->type() == MAST::STRUCTURAL_COMPLIANCE) {
            
            libMesh::NonlinearImplicitSystem& sys =
            dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
            
            MAST::RealOutputFunction& output =
            dynamic_cast<MAST::RealOutputFunction&>(*it->second);
            
            for ( unsigned int i=0; i<params.size(); i++) {
                
                const MAST::FunctionBase* f =
                _discipline->get_parameter(&(params[i].get()));
                
                if (!if_total_sensitivity)
                    _calculate_compliance(X,
                                          sys,
                                          output,
                                          f);
                else {
                    
                    localized_solution_sensitivity.reset
                    (_build_localized_vector(sys,
                                             sys.get_sensitivity_solution(i)).release());
                    
                    
                    _calculate_compliance(X,
                                          sys,
                                          output,
                                          f,
                                          localized_solution_sensitivity.get());
                }
            }
        }
    
    // now call the parent's method for calculation of the other outputs
    MAST::NonlinearImplicitAssembly::
    calculate_output_sensitivity(params,
                                 if_total_sensitivity,
                                 X);
    
}






void
MAST::StructuralNonlinearAssembly::
_calculate_compliance (const libMesh::NumericVector<Real>& X,
                       libMesh::NonlinearImplicitSystem& S,
                       MAST::RealOutputFunction& output,
                       const MAST::FunctionBase* f,
                       const libMesh::NumericVector<Real>* dX) {


    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &(nonlin_sys));
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol, dsol;
    RealMatrixX mat;
    Real        val;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_solution_sens;
    localized_solution.reset(_build_localized_vector(nonlin_sys,
                                                     X).release());
    if (f)
        localized_solution_sens.reset(_build_localized_vector(nonlin_sys,
                                                              *dX).release());

    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        physics_elem.reset(_build_elem(*elem).release());
        
        MAST::StructuralElementBase& p_elem =
        dynamic_cast<MAST::StructuralElementBase&>(*physics_elem);
        
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero (ndofs);
        dsol.setZero(ndofs);
        vec.setZero (ndofs);
        mat.setZero (ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);

        physics_elem->set_solution(sol, false);  // primal solution

        if (f) {
            
            for (unsigned int i=0; i<dof_indices.size(); i++)
                dsol(i) = (*localized_solution_sens)(dof_indices[i]);
            
            physics_elem->set_solution(dsol, true);  // sensitivity solution
            physics_elem->sensitivity_param = f;
        }
        
        // set the incompatible mode solution if required by the
        // element
        if (p_elem.if_incompatible_modes()) {
            // check if the vector exists in the map
            if (!_incompatible_sol.count(elem))
                _incompatible_sol[elem] = RealVectorX::Zero(p_elem.incompatible_mode_size());
            p_elem.set_incompatible_mode_solution(_incompatible_sol[elem]);
        }
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);

        // if no parameter has been specified, then only evaluate
        // the output. Otherwise, calculate the sensitivity
        if (!f) {
            
            // perform the element level calculations
            _elem_calculations(*physics_elem,
                               true,
                               vec, mat);
            
            
            output.add_value(sol.dot(mat * sol));
        }
        else {
            
            // perform the element level calculations
            _elem_calculations(*physics_elem,
                               true,
                               vec, mat);
            val = sol.dot(mat * dsol) + dsol.dot(mat *  sol);
            output.add_sensitivity(f, val);
            
            _elem_sensitivity_calculations(*physics_elem, true, vec, mat);
            output.add_sensitivity(f, sol.dot(mat * sol));
        }
        
        physics_elem->detach_active_solution_function();
        
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
}




std::auto_ptr<MAST::ElementBase>
MAST::StructuralNonlinearAssembly::_build_elem(const libMesh::Elem& elem) {
    
    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    MAST::ElementBase* rval =
    MAST::build_structural_element(*_system, elem, p).release();
    
    return std::auto_ptr<MAST::ElementBase>(rval);
}




void
MAST::StructuralNonlinearAssembly::
_elem_calculations(MAST::ElementBase& elem,
                   bool if_jac,
                   RealVectorX& vec,
                   RealMatrixX& mat) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    vec.setZero();
    mat.setZero();
    RealMatrixX
    dummy = RealMatrixX::Zero(mat.rows(), mat.cols());
    
    e.internal_residual(if_jac, vec, mat);
    e.side_external_residual(if_jac,
                             vec,
                             dummy,
                             mat,
                             _discipline->side_loads());
    e.volume_external_residual(if_jac,
                               vec,
                               dummy,
                               mat,
                               _discipline->volume_loads());
}



void
MAST::StructuralNonlinearAssembly::
_elem_linearized_jacobian_solution_product(MAST::ElementBase& elem,
                                           RealVectorX& vec) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    vec.setZero();

    libmesh_error(); // to be implemented
}




void
MAST::StructuralNonlinearAssembly::
_elem_sensitivity_calculations(MAST::ElementBase& elem,
                               bool if_jac,
                               RealVectorX& vec,
                               RealMatrixX& mat) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    vec.setZero();
    mat.setZero();
    RealMatrixX
    dummy = RealMatrixX::Zero(vec.size(), vec.size());
    
    e.internal_residual_sensitivity(if_jac, vec, mat);
    e.side_external_residual_sensitivity(if_jac,
                                         vec,
                                         dummy,
                                         mat,
                                         _discipline->side_loads());
    e.volume_external_residual_sensitivity(if_jac,
                                           vec,
                                           dummy,
                                           mat,
                                           _discipline->volume_loads());
}




void
MAST::StructuralNonlinearAssembly::
attach_discipline_and_system(MAST::PhysicsDisciplineBase& discipline,
                             MAST::SystemInitialization& system) {
    
    // call the parent's method firts
    MAST::NonlinearImplicitAssembly::attach_discipline_and_system(discipline,
                                                                  system);
    
    // get the nonlinear solver SNES object from System and
    // add a monitor to it so that it can be used to update the
    // incompatible mode solution after each update
    
    libMesh::PetscNonlinearSolver<Real> &petsc_nonlinear_solver =
    *(dynamic_cast<libMesh::PetscNonlinearSolver<Real>*>
      (dynamic_cast<libMesh::NonlinearImplicitSystem&>(system.system()).nonlinear_solver.get()));


    // initialize the solver before getting the snes object
    if (libMesh::on_command_line("--solver_system_names"))
        petsc_nonlinear_solver.init((system.system().name()+"_").c_str());
    
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
MAST::StructuralNonlinearAssembly::
clear_discipline_and_system( ) {
    
    // next, remove the monitor function from the snes object
    libMesh::PetscNonlinearSolver<Real> &petsc_nonlinear_solver =
    *(dynamic_cast<libMesh::PetscNonlinearSolver<Real>*>
      (dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system()).nonlinear_solver.get()));
    
    // get the SNES object
    SNES snes = petsc_nonlinear_solver.snes();
    
    PetscErrorCode ierr =
    SNESMonitorCancel(snes);
    libmesh_assert(!ierr);
    
    // call the parent's method firts
    MAST::NonlinearImplicitAssembly::clear_discipline_and_system();
    
    
}



