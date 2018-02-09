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
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_assembly.h"
#include "property_cards/element_property_card_1D.h"
#include "base/physics_discipline_base.h"
#include "mesh/local_elem_fe.h"




MAST::StructuralNonlinearAssemblyElemOperations::StructuralNonlinearAssemblyElemOperations():
MAST::NonlinearImplicitAssemblyElemOperations(),
_incompatible_sol_assembly(nullptr) {
    
}



MAST::StructuralNonlinearAssemblyElemOperations::
~StructuralNonlinearAssemblyElemOperations() {
    
}



void
MAST::StructuralNonlinearAssemblyElemOperations::set_elem_solution(const RealVectorX& sol) {
    
    unsigned int
    n = (unsigned int)sol.size();
    
    RealVectorX
    zero = RealVectorX::Zero(n);
    
    _physics_elem->set_solution    (sol);
    _physics_elem->set_velocity    (zero); // set to zero vector for a quasi-steady analysis
    _physics_elem->set_acceleration(zero); // set to zero vector for a quasi-steady analysis
    

    if (_incompatible_sol_assembly) {
        
        // set the incompatible mode solution if required by the
        // element
        
        MAST::StructuralElementBase& s_elem =
        dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        
        if (s_elem.if_incompatible_modes())
            _incompatible_sol_assembly->set_elem_incompatible_sol(s_elem);
    }
}



void
MAST::StructuralNonlinearAssemblyElemOperations::
init(const libMesh::Elem& elem) {
    
    libmesh_assert(!_physics_elem);
    libmesh_assert(_system);
    libmesh_assert(_assembly);

    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>
    (_discipline->get_property_card(elem));
    
    _physics_elem =
    MAST::build_structural_element(*_system, *_assembly, elem, p).release();
}




void
MAST::StructuralNonlinearAssemblyElemOperations::
elem_calculations(bool if_jac,
                  RealVectorX& vec,
                  RealMatrixX& mat) {
    
    libmesh_assert(_physics_elem);
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    
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
MAST::StructuralNonlinearAssemblyElemOperations::
elem_linearized_jacobian_solution_product(RealVectorX& vec) {
    
    libmesh_assert(_physics_elem);

    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    
    vec.setZero();
    RealMatrixX
    mat = RealMatrixX::Zero(vec.size(), vec.size());
    
    e.linearized_internal_residual(false, vec, mat);
    e.linearized_side_external_residual(false,
                                        vec,
                                        mat,
                                        mat,
                                        _discipline->side_loads());
    e.linearized_volume_external_residual(false,
                                          vec,
                                          mat,
                                          mat,
                                          _discipline->volume_loads());
}




void
MAST::StructuralNonlinearAssemblyElemOperations::
elem_sensitivity_calculations(RealVectorX& vec) {
    
    libmesh_assert(_physics_elem);

    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    
    vec.setZero();
    RealMatrixX
    dummy = RealMatrixX::Zero(vec.size(), vec.size());
    
    e.internal_residual_sensitivity(false, vec, dummy);
    e.side_external_residual_sensitivity(false,
                                         vec,
                                         dummy,
                                         dummy,
                                         _discipline->side_loads());
    e.volume_external_residual_sensitivity(false,
                                           vec,
                                           dummy,
                                           dummy,
                                           _discipline->volume_loads());
}



void
MAST::StructuralNonlinearAssemblyElemOperations::
elem_second_derivative_dot_solution_assembly(RealMatrixX& m) {
    
    libmesh_assert(_physics_elem);
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    m.setZero();
    
    e.internal_residual_jac_dot_state_sensitivity(m);
}




//void
//MAST::StructuralNonlinearAssemblyElemOperations::
//calculate_outputs(const libMesh::NumericVector<Real>& X) {
//
//
//    // look through the outputs for structural compliance. If present,
//    // then evaluate it separately
//    MAST::VolumeOutputMapType::iterator
//    it    = _discipline->volume_output().begin(),
//    end   = _discipline->volume_output().end();
//
//
//    for ( ; it != end; it++)
//        if (it->second->type() == MAST::STRUCTURAL_COMPLIANCE) {
//
//            MAST::NonlinearSystem& sys = _system->system();
//
//            MAST::RealOutputFunction& output =
//            dynamic_cast<MAST::RealOutputFunction&>(*it->second);
//
//            _calculate_compliance(X, sys, output);
//        }
//
//    // now call the parent's method for calculation of the other outputs
//    MAST::NonlinearImplicitAssembly::calculate_outputs(X);
//
//}
//
//
//
//
//
//
//void
//MAST::StructuralNonlinearAssemblyElemOperations::
//calculate_output_sensitivity(libMesh::ParameterVector& params,
//                             const bool if_total_sensitivity,
//                             const libMesh::NumericVector<Real>& X) {
//
//
//    // look through the outputs for structural compliance. If present,
//    // then evaluate it separately
//    MAST::VolumeOutputMapType::iterator
//    it    = _discipline->volume_output().begin(),
//    end   = _discipline->volume_output().end();
//
//    std::unique_ptr<libMesh::NumericVector<Real> >
//    localized_solution_sensitivity;
//
//    for ( ; it != end; it++)
//        if (it->second->type() == MAST::STRUCTURAL_COMPLIANCE) {
//
//            MAST::NonlinearSystem& sys = _system->system();
//
//            MAST::RealOutputFunction& output =
//            dynamic_cast<MAST::RealOutputFunction&>(*it->second);
//
//            for ( unsigned int i=0; i<params.size(); i++) {
//
//                const MAST::FunctionBase* f =
//                _discipline->get_parameter(&(params[i].get()));
//
//                if (!if_total_sensitivity)
//                    _calculate_compliance(X,
//                                          sys,
//                                          output,
//                                          f);
//                else {
//
//                    localized_solution_sensitivity.reset
//                    (build_localized_vector(sys,
//                                             sys.get_sensitivity_solution(i)).release());
//
//
//                    _calculate_compliance(X,
//                                          sys,
//                                          output,
//                                          f,
//                                          localized_solution_sensitivity.get());
//                }
//            }
//        }
//
//    // now call the parent's method for calculation of the other outputs
//    MAST::NonlinearImplicitAssembly::
//    calculate_output_sensitivity(params,
//                                 if_total_sensitivity,
//                                 X);
//
//}






//void
//MAST::StructuralNonlinearAssemblyElemOperations::
//_calculate_compliance (const libMesh::NumericVector<Real>& X,
//                       libMesh::NonlinearImplicitSystem& S,
//                       MAST::RealOutputFunction& output,
//                       const MAST::FunctionBase* f,
//                       const libMesh::NumericVector<Real>* dX) {
//
//
//    MAST::NonlinearSystem& nonlin_sys = _system->system();
//
//    // make sure that the system for which this object was created,
//    // and the system passed through the function call are the same
//    libmesh_assert_equal_to(&S, &(nonlin_sys));
//
//    // iterate over each element, initialize it and get the relevant
//    // analysis quantities
//    RealVectorX vec, sol, dsol;
//    RealMatrixX mat;
//    Real        val;
//
//    std::vector<libMesh::dof_id_type> dof_indices;
//    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
//    
//
//    std::unique_ptr<libMesh::NumericVector<Real> >
//    localized_solution,
//    localized_solution_sens;
//    localized_solution.reset(build_localized_vector(nonlin_sys,
//                                                     X).release());
//    if (f)
//        localized_solution_sens.reset(build_localized_vector(nonlin_sys,
//                                                              *dX).release());
//
//
//    // if a solution function is attached, initialize it
//    if (_sol_function)
//        _sol_function->init( X);
//
//
//    libMesh::MeshBase::const_element_iterator       el     =
//    nonlin_sys.get_mesh().active_local_elements_begin();
//    const libMesh::MeshBase::const_element_iterator end_el =
//    nonlin_sys.get_mesh().active_local_elements_end();
//
//
//    for ( ; el != end_el; ++el) {
//
//        const libMesh::Elem* elem = *el;
//
//        dof_map.dof_indices (elem, dof_indices);
//
//        _build_elem(*elem).release());
//
//        MAST::StructuralElementBase& p_elem =
//        dynamic_cast<MAST::StructuralElementBase&>(*physics_elem);
//
//
//        // get the solution
//        unsigned int ndofs = (unsigned int)dof_indices.size();
//        sol.setZero (ndofs);
//        dsol.setZero(ndofs);
//        vec.setZero (ndofs);
//        mat.setZero (ndofs, ndofs);
//
//        for (unsigned int i=0; i<dof_indices.size(); i++)
//            sol(i) = (*localized_solution)(dof_indices[i]);
//
//        physics_elem->set_solution(sol, false);  // primal solution
//
//        if (f) {
//
//            for (unsigned int i=0; i<dof_indices.size(); i++)
//                dsol(i) = (*localized_solution_sens)(dof_indices[i]);
//
//            physics_elem->set_solution(dsol, true);  // sensitivity solution
//            physics_elem->sensitivity_param = f;
//        }
//
//        // set the incompatible mode solution if required by the
//        // element
//        if (p_elem.if_incompatible_modes()) {
//            // check if the vector exists in the map
//            if (!_incompatible_sol.count(elem))
//                _incompatible_sol[elem] = RealVectorX::Zero(p_elem.incompatible_mode_size());
//            p_elem.set_incompatible_mode_solution(_incompatible_sol[elem]);
//        }
//
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
//
//        // if no parameter has been specified, then only evaluate
//        // the output. Otherwise, calculate the sensitivity
//        if (!f) {
//
//            // perform the element level calculations
//            _elem_calculations(*physics_elem,
//                               true,
//                               vec, mat);
//
//
//            output.add_value(sol.dot(mat * sol));
//        }
//        else {
//
//            // perform the element level calculations
//            _elem_calculations(*physics_elem,
//                               true,
//                               vec, mat);
//            val = sol.dot(mat * dsol) + dsol.dot(mat *  sol);
//            output.add_sensitivity(f, val);
//
//            _elem_sensitivity_calculations(*physics_elem, true, vec, mat);
//            output.add_sensitivity(f, sol.dot(mat * sol));
//        }
//
//        physics_elem->detach_active_solution_function();
//
//    }
//
//
//    // if a solution function is attached, clear it
//    if (_sol_function)
//        _sol_function->clear();
//}



void
MAST::StructuralNonlinearAssemblyElemOperations::
set_local_fe_data(MAST::LocalElemFE& fe,
                  const libMesh::Elem& e) const {
     
    if (e.dim() == 1) {
        
        const MAST::ElementPropertyCard1D&
        p_card = dynamic_cast<const MAST::ElementPropertyCard1D&>
        (_discipline->get_property_card(e));
        
        fe.set_1d_y_vector(p_card.y_vector());
    }
}







