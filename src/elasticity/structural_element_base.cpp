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
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_element_1d.h"
#include "elasticity/structural_element_2d.h"
#include "elasticity/solid_element_3d.h"
#include "base/system_initialization.h"
#include "base/boundary_condition_base.h"
#include "property_cards/element_property_card_1D.h"
#include "mesh/local_elem_base.h"
#include "mesh/local_1d_elem.h"
#include "mesh/local_2d_elem.h"
#include "mesh/local_3d_elem.h"
#include "numerics/fem_operator_matrix.h"
#include "numerics/utility.h"
#include "elasticity/stress_output_base.h"



MAST::StructuralElementBase::StructuralElementBase(MAST::SystemInitialization& sys,
                                                   const libMesh::Elem& elem,
                                                   const MAST::ElementPropertyCardBase& p):
MAST::ElementBase(sys, elem),
follower_forces(false),
_property(p),
_incompatible_sol(NULL) {
    
    MAST::LocalElemBase* rval = NULL;
    
    switch (elem.dim()) {
        case 1: {
            const MAST::ElementPropertyCard1D& p_1d =
            dynamic_cast<const MAST::ElementPropertyCard1D&>(p);
            rval = new MAST::Local1DElem(elem, p_1d.y_vector());
        }
            break;
            
        case 2:
            rval = new MAST::Local2DElem(elem);
            break;
            
        case 3:
            rval = new MAST::Local3DElem(elem);
            break;
            
        default:
            // should not get here.
            libmesh_error();
            break;
    }
    
    _local_elem.reset(rval);
    
}



MAST::StructuralElementBase::~StructuralElementBase() {
    
}



void
MAST::StructuralElementBase::set_solution(const RealVectorX& vec,
                                          bool if_sens) {
    
    // convert the vector to the local element coordinate system
    if (!if_sens) {
        if (_elem.dim() == 3)
            _local_sol = vec;
        else {
            _local_sol = RealVectorX::Zero(vec.size());
            this->transform_vector_to_local_system(vec, _local_sol);
        }
    }
    else {
        
        // set the element solution sensitivity.
        if (_elem.dim() == 3)
            _local_sol_sens = vec;
        else {
            _local_sol_sens = RealVectorX::Zero(vec.size());
            this->transform_vector_to_local_system(vec, _local_sol_sens);
        }
    }
    
    MAST::ElementBase::set_solution(vec, if_sens);
}



void
MAST::StructuralElementBase::set_velocity(const RealVectorX& vec,
                                          bool if_sens) {
    
    if (!if_sens) {
        if (_elem.dim() == 3)
            _local_vel = vec;
        else {
            _local_vel = RealVectorX::Zero(vec.size());
            this->transform_vector_to_local_system(vec, _local_vel);
        }
    }
    else {
        
        if (_elem.dim() == 3)
            _local_vel_sens = vec;
        else {
            _local_vel_sens = RealVectorX::Zero(vec.size());
            this->transform_vector_to_local_system(vec, _local_vel_sens);
        }
    }

    MAST::ElementBase::set_velocity(vec, if_sens);
}



void
MAST::StructuralElementBase::set_acceleration(const RealVectorX& vec,
                                              bool if_sens) {
    
    if (!if_sens) {
        if (_elem.dim() == 3)
            _local_accel = vec;
        else {
            _local_accel = RealVectorX::Zero(vec.size());
            this->transform_vector_to_local_system(vec, _local_accel);
        }
    }
    else {

        if (_elem.dim() == 3)
            _local_accel_sens = vec;
        else {
            _local_accel_sens = RealVectorX::Zero(vec.size());
            this->transform_vector_to_local_system(vec, _local_accel_sens);
        }
    }
    
    MAST::ElementBase::set_acceleration(vec, if_sens);
}



//void
//MAST::StructuralElementBase::set_base_solution(const RealVectorX& vec,
//                                               bool if_sens) {
//    
//    if (!if_sens) {
//        if (_elem.dim() == 3)
//            _local_base_sol = vec;
//        else {
//            _local_base_sol = RealVectorX::Zero(vec.size());
//            this->transform_vector_to_local_system(vec, _local_base_sol);
//        }
//    }
//    else
//        libmesh_error();
//    
//    MAST::ElementBase::set_base_solution(vec, if_sens);
//}



bool
MAST::StructuralElementBase::inertial_residual (bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac_xddot,
                                                RealMatrixX& jac_xdot,
                                                RealMatrixX& jac) {
    
    const std::vector<Real>& JxW               = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz     = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = _fe->get_phi();
    
    const unsigned int
    n_phi    = (unsigned int)phi.size(),
    n_vars   = _system.system().n_vars(),
    n1       =6,
    n2       =6*n_phi;
    
    RealMatrixX
    material_mat,
    mat1_n1n2     = RealMatrixX::Zero(n1, n2),
    mat2_n2n2     = RealMatrixX::Zero(n2, n2),
    local_jac     = RealMatrixX::Zero(n2, n2);
    RealVectorX
    phi_vec    = RealVectorX::Zero(n_phi),
    vec1_n1    = RealVectorX::Zero(n1),
    vec2_n2    = RealVectorX::Zero(n2),
    local_f    = RealVectorX::Zero(n2);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
    mat_inertia  = _property.inertia_matrix(*this);
    
    libMesh::Point p;
    MAST::FEMOperatorMatrix Bmat;
    
    if (_property.if_diagonal_mass_matrix()) {
        
        // as an approximation, get matrix at the first quadrature point
        _local_elem->global_coordinates_location(xyz[0], p);
        
        (*mat_inertia)(p, _time, material_mat);
        
        Real vol = 0.;
        const unsigned int nshp = _fe->n_shape_functions();
        for (unsigned int i=0; i<JxW.size(); i++)
            vol += JxW[i];
        vol /= (1.* nshp);
        for (unsigned int i_var=0; i_var<6; i_var++)
            for (unsigned int i=0; i<nshp; i++)
                local_jac(i_var*nshp+i, i_var*nshp+i) =
                vol*material_mat(i_var, i_var);
        
        local_f =  local_jac * _local_accel;
    }
    else {
        libMesh::Point p;
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            
            _local_elem->global_coordinates_location(xyz[0], p);
            
            (*mat_inertia)(p, _time, material_mat);
            
            // now set the shape function values
            for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
                phi_vec(i_nd) = phi[i_nd][qp];
            
            Bmat.reinit(n_vars, phi_vec);
            
            Bmat.left_multiply(mat1_n1n2, material_mat);
            
            vec1_n1 = mat1_n1n2 * _local_accel;
            Bmat.vector_mult_transpose(vec2_n2, vec1_n1);
            
            local_f += JxW[qp] * vec2_n2;
            
            if (request_jacobian) {
                
                Bmat.right_multiply_transpose(mat2_n2n2,
                                              mat1_n1n2);
                local_jac += JxW[qp]*mat2_n2n2;
            }
        }
    }
    
    // now transform to the global coorodinate system
    if (_elem.dim() < 3) {
        transform_vector_to_global_system(local_f, vec2_n2);
        f += vec2_n2;
        
        if (request_jacobian) {
            transform_matrix_to_global_system(local_jac, mat2_n2n2);
            jac_xddot += mat2_n2n2;
        }
    }
    else {
        
        f += local_f;
        if (request_jacobian)
            jac_xddot += local_jac;
    }
    
    return request_jacobian;
}







bool
MAST::StructuralElementBase::inertial_residual_sensitivity (bool request_jacobian,
                                                            RealVectorX& f,
                                                            RealMatrixX& jac_xddot,
                                                            RealMatrixX& jac_xdot,
                                                            RealMatrixX& jac) {
    
    const std::vector<Real>& JxW               = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz     = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = _fe->get_phi();
    
    const unsigned int
    n_phi    = (unsigned int)phi.size(),
    n_vars   = _system.system().n_vars(),
    n1       =6,
    n2       =6*n_phi;
    
    RealMatrixX
    material_mat,
    mat1_n1n2     = RealMatrixX::Zero(n1, n2),
    mat2_n2n2     = RealMatrixX::Zero(n2, n2),
    local_jac     = RealMatrixX::Zero(n2, n2);
    RealVectorX
    phi_vec    = RealVectorX::Zero(n_phi),
    vec1_n1    = RealVectorX::Zero(n1),
    vec2_n2    = RealVectorX::Zero(n2),
    local_f    = RealVectorX::Zero(n2);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
    mat_inertia  = _property.inertia_matrix(*this);
    
    libMesh::Point p;
    MAST::FEMOperatorMatrix Bmat;
    
    if (_property.if_diagonal_mass_matrix()) {
        
        // as an approximation, get matrix at the first quadrature point
        _local_elem->global_coordinates_location(xyz[0], p);
        
        mat_inertia->derivative(MAST::PARTIAL_DERIVATIVE,
                                *this->sensitivity_param,
                                p, _time, material_mat);
        
        Real vol = 0.;
        const unsigned int nshp = _fe->n_shape_functions();
        for (unsigned int i=0; i<JxW.size(); i++)
            vol += JxW[i];
        vol /= (1.* nshp);
        for (unsigned int i_var=0; i_var<6; i_var++)
            for (unsigned int i=0; i<nshp; i++)
                local_jac(i_var*nshp+i, i_var*nshp+i) =
                vol*material_mat(i_var, i_var);
        
        local_f =  local_jac * _local_accel;
    }
    else {
        libMesh::Point p;
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            
            _local_elem->global_coordinates_location(xyz[0], p);
            
            mat_inertia->derivative(MAST::PARTIAL_DERIVATIVE,
                                    *this->sensitivity_param,
                                    p, _time, material_mat);
            
            // now set the shape function values
            for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
                phi_vec(i_nd) = phi[i_nd][qp];
            
            Bmat.reinit(n_vars, phi_vec);
            
            Bmat.left_multiply(mat1_n1n2, material_mat);
            
            vec1_n1 = mat1_n1n2 * _local_accel;
            Bmat.vector_mult_transpose(vec2_n2, vec1_n1);
            
            local_f += JxW[qp] * vec2_n2;
            
            if (request_jacobian) {
                
                Bmat.right_multiply_transpose(mat2_n2n2,
                                              mat1_n1n2);
                local_jac += JxW[qp]*mat2_n2n2;
            }
        }
    }
    
    // now transform to the global coorodinate system
    if (_elem.dim() < 3) {
        transform_vector_to_global_system(local_f, vec2_n2);
        f += vec2_n2;
        
        if (request_jacobian) {
            transform_matrix_to_global_system(local_jac, mat2_n2n2);
            jac_xddot += mat2_n2n2;
        }
    }
    else {
        
        f += local_f;
        if (request_jacobian)
            jac_xddot += local_jac;
    }
    
    return request_jacobian;
}





template <typename ValType>
bool
MAST::StructuralElementBase::
side_external_residual(bool request_jacobian,
                       RealVectorX& f,
                       RealMatrixX& jac_xdot,
                       RealMatrixX& jac,
                       std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    typedef std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    const libMesh::BoundaryInfo& binfo = *_system.system().get_mesh().boundary_info;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    for (unsigned short int n=0; n<_elem.n_sides(); n++) {
        
        // if no boundary ids have been specified for the side, then
        // move to the next side.
        if (!binfo.n_boundary_ids(&_elem, n))
            continue;
        
        // check to see if any of the specified boundary ids has a boundary
        // condition associated with them
        std::vector<libMesh::boundary_id_type> bc_ids = binfo.boundary_ids(&_elem, n);
        std::vector<libMesh::boundary_id_type>::const_iterator bc_it = bc_ids.begin();
        
        for ( ; bc_it != bc_ids.end(); bc_it++) {
            
            // find the loads on this boundary and evaluate the f and jac
            it = bc.equal_range(*bc_it);
            
            for ( ; it.first != it.second; it.first++) {
                
                // apply all the types of loading
                switch (it.first->second->type()) {
                    case MAST::SURFACE_PRESSURE:
                        calculate_jac = (calculate_jac ||
                                         surface_pressure_residual(request_jacobian,
                                                                   f, jac,
                                                                   n,
                                                                   *it.first->second));
                        break;

                        
                    case MAST::PISTON_THEORY:
                        calculate_jac = (calculate_jac ||
                                         piston_theory_residual(request_jacobian,
                                                                f,
                                                                jac_xdot,
                                                                jac,
                                                                n,
                                                                *it.first->second));
                        break;

                        
                    case MAST::SMALL_DISTURBANCE_MOTION:
                        calculate_jac = (calculate_jac ||
                                         small_disturbance_surface_pressure_residual<ValType>
                                         (request_jacobian,
                                          f, jac,
                                          n,
                                          *it.first->second));
                        break;
                        
                        
                    case MAST::DIRICHLET:
                        // nothing to be done here
                        break;
                        
                    default:
                        // not implemented yet
                        libmesh_error();
                        break;
                }
            }
        }
    }
    return (request_jacobian && calculate_jac);
}


template <typename ValType>
bool
MAST::StructuralElementBase::
volume_external_residual (bool request_jacobian,
                          RealVectorX& f,
                          RealMatrixX& jac_xdot,
                          RealMatrixX& jac,
                          std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    // iterate over the boundary ids given in the provided force map
    std::pair<std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator,
    std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    libMesh::subdomain_id_type sid = _elem.subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::SURFACE_PRESSURE:
                calculate_jac = (calculate_jac ||
                                 surface_pressure_residual(request_jacobian,
                                                           f, jac,
                                                           *it.first->second));
                break;

            case MAST::PISTON_THEORY:
                calculate_jac = (calculate_jac ||
                                 piston_theory_residual(request_jacobian,
                                                        f,
                                                        jac_xdot,
                                                        jac,
                                                        *it.first->second));
                break;

            case MAST::TEMPERATURE:
                calculate_jac = (calculate_jac ||
                                 thermal_residual(request_jacobian,
                                                  f, jac,
                                                  *it.first->second));
                break;
                
            case MAST::SMALL_DISTURBANCE_MOTION:
                calculate_jac = (calculate_jac ||
                                 small_disturbance_surface_pressure_residual<ValType>
                                 (request_jacobian,
                                  f, jac,
                                  *it.first->second));
                break;
                
            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    
    return (request_jacobian && calculate_jac);
}






template <typename ValType>
bool
MAST::StructuralElementBase::
side_external_residual_sensitivity(bool request_jacobian,
                                   RealVectorX& f,
                                   RealMatrixX& jac_xdot,
                                   RealMatrixX& jac,
                                   std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    typedef std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    const libMesh::BoundaryInfo& binfo = *_system.system().get_mesh().boundary_info;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    for (unsigned short int n=0; n<_elem.n_sides(); n++) {
        
        // if no boundary ids have been specified for the side, then
        // move to the next side.
        if (!binfo.n_boundary_ids(&_elem, n))
            continue;
        
        // check to see if any of the specified boundary ids has a boundary
        // condition associated with them
        std::vector<libMesh::boundary_id_type> bc_ids = binfo.boundary_ids(&_elem, n);
        std::vector<libMesh::boundary_id_type>::const_iterator bc_it = bc_ids.begin();
        
        for ( ; bc_it != bc_ids.end(); bc_it++) {
            
            // find the loads on this boundary and evaluate the f and jac
            it = bc.equal_range(*bc_it);
            
            for ( ; it.first != it.second; it.first++) {
                
                // apply all the types of loading
                switch (it.first->second->type()) {
                    case MAST::SURFACE_PRESSURE:
                        calculate_jac = (calculate_jac ||
                                         surface_pressure_residual_sensitivity(request_jacobian,
                                                                               f, jac,
                                                                               n,
                                                                               *it.first->second));
                        break;
                        
                        
                    case MAST::PISTON_THEORY:
                        calculate_jac = (calculate_jac ||
                                         piston_theory_residual_sensitivity(request_jacobian,
                                                                            f,
                                                                            jac_xdot,
                                                                            jac,
                                                                            n,
                                                                            *it.first->second));
                        break;
                        
                        
                    case MAST::SMALL_DISTURBANCE_MOTION:
                        /*
                        calculate_jac = (calculate_jac ||
                                         small_disturbance_surface_pressure_residual_sensitivity<ValType>
                                         (request_jacobian,
                                          f, jac,
                                          n,
                                          *it.first->second));
                         */
                        libmesh_error(); // to be implemented
                        break;
                        
                        
                    case MAST::DIRICHLET:
                        // nothing to be done here
                        break;
                        
                    default:
                        // not implemented yet
                        libmesh_error();
                        break;
                }
            }
        }
    }
    return (request_jacobian && calculate_jac);
}


template <typename ValType>
bool
MAST::StructuralElementBase::
volume_external_residual_sensitivity (bool request_jacobian,
                                      RealVectorX& f,
                                      RealMatrixX& jac_xdot,
                                      RealMatrixX& jac,
                                      std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    // iterate over the boundary ids given in the provided force map
    std::pair<std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator,
    std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    libMesh::subdomain_id_type sid = _elem.subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::SURFACE_PRESSURE:
                calculate_jac = (calculate_jac ||
                                 surface_pressure_residual_sensitivity(request_jacobian,
                                                                       f, jac,
                                                                       *it.first->second));
                break;
                
            case MAST::PISTON_THEORY:
                calculate_jac = (calculate_jac ||
                                 piston_theory_residual_sensitivity(request_jacobian,
                                                                    f,
                                                                    jac_xdot,
                                                                    jac,
                                                                    *it.first->second));
                break;
                
            case MAST::TEMPERATURE:
                calculate_jac = (calculate_jac ||
                                 thermal_residual_sensitivity(request_jacobian,
                                                              f, jac,
                                                              *it.first->second));
                break;
                
            case MAST::SMALL_DISTURBANCE_MOTION:
                /*
                calculate_jac = (calculate_jac ||
                                 small_disturbance_surface_pressure_residual_sensitivity<ValType>
                                 (request_jacobian,
                                  f, jac,
                                  *it.first->second));
                 */
                libmesh_error(); // to be implemented
                break;
                
            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    
    return (request_jacobian && calculate_jac);
}





bool
MAST::StructuralElementBase::volume_output_quantity
(bool request_derivative,
 bool request_sensitivity,
 std::multimap<libMesh::subdomain_id_type, MAST::OutputFunctionBase*>& output) {
    
    
    
    // iterate over the boundary ids given in the provided force map
    std::pair<std::multimap<libMesh::subdomain_id_type, MAST::OutputFunctionBase*>::const_iterator,
    std::multimap<libMesh::subdomain_id_type, MAST::OutputFunctionBase*>::const_iterator> it;
    
    // for each subdomain id, check if the element has this domain id
    bool calculate_derivative = false;
    
    libMesh::subdomain_id_type sid = _elem.subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =  output.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::STRAIN_STRESS_TENSOR: {
                // evaluate the stress only if the object needs it
                MAST::StressStrainOutputBase& stress =
                dynamic_cast<MAST::StressStrainOutputBase&>(*(it.first->second));
                if (stress.evaluate_for_element(_elem)) {
                    // look through the discipline to see if there are
                    // any thermal stress conditions defined for this
                    // element.
                    calculate_derivative = (calculate_derivative ||
                                            calculate_stress(request_derivative,
                                                             request_sensitivity,
                                                             *it.first->second));
                }
            }
                break;
                
            case MAST::STRUCTURAL_COMPLIANCE:
                // nothing to be done here since it is handled separately
                // in the structural implicit assembly
                break;
                
            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    
    return ((request_derivative || request_sensitivity) && calculate_derivative);

    return false;
}






bool
MAST::StructuralElementBase::
surface_pressure_residual(bool request_jacobian,
                          RealVectorX &f,
                          RealMatrixX &jac,
                          MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    
    const std::vector<Real> &JxW                = _fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint   = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi  = _fe->get_phi();
    const unsigned int
    n_phi = (unsigned int)phi.size(),
    n1    =3,
    n2    =6*n_phi;
    
    
    // normal for face integration
    libMesh::Point normal;
    // direction of pressure assumed to be normal (along local z-axis)
    // to the element face for 2D and along local y-axis for 1D element.
    normal(_elem.dim()) = -1.;
    
    
    // get the function from this boundary condition
    MAST::FieldFunction<Real>& func =
    bc.get<MAST::FieldFunction<Real> >("pressure");
    
    Real press;
    FEMOperatorMatrix Bmat;
    libMesh::Point pt;
    
    RealVectorX
    phi_vec  = RealVectorX::Zero(n_phi),
    force    = RealVectorX::Zero(2*n1),
    local_f  = RealVectorX::Zero(n2),
    vec_n2   = RealVectorX::Zero(n2);
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++)
    {
        
        _local_elem->global_coordinates_location(qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure value
        func(pt, _time, press);
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = press * normal(i_dim);
        
        Bmat.vector_mult_transpose(vec_n2, force);
        
        local_f += JxW[qp] * vec_n2;
    }
    
    // now transform to the global system and add
    transform_vector_to_global_system(local_f, vec_n2);
    f -= vec_n2;
    
    
    return (request_jacobian);
}






bool
MAST::StructuralElementBase::
surface_pressure_residual_sensitivity(bool request_jacobian,
                                      RealVectorX &f,
                                      RealMatrixX &jac,
                                      MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    
    const std::vector<Real> &JxW                = _fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint   = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi  = _fe->get_phi();
    const unsigned int
    n_phi = (unsigned int)phi.size(),
    n1    =3,
    n2    =6*n_phi;
    
    
    // normal for face integration
    libMesh::Point normal;
    // direction of pressure assumed to be normal (along local z-axis)
    // to the element face for 2D and along local y-axis for 1D element.
    normal(_elem.dim()) = -1.;
    
    
    // get the function from this boundary condition
    MAST::FieldFunction<Real>& func =
    bc.get<MAST::FieldFunction<Real> >("pressure");
    
    Real press;
    FEMOperatorMatrix Bmat;
    libMesh::Point pt;
    
    RealVectorX
    phi_vec  = RealVectorX::Zero(n_phi),
    force    = RealVectorX::Zero(2*n1),
    local_f  = RealVectorX::Zero(n2),
    vec_n2   = RealVectorX::Zero(n2);
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++)
    {
        
        _local_elem->global_coordinates_location(qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure value
        func.derivative(MAST::TOTAL_DERIVATIVE,
                        *sensitivity_param,
                        pt,
                        _time,
                        press);
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = press * normal(i_dim);
        
        Bmat.vector_mult_transpose(vec_n2, force);
        
        local_f += JxW[qp] * vec_n2;
    }
    
    // now transform to the global system and add
    transform_vector_to_global_system(local_f, vec_n2);
    f -= vec_n2;
    
    
    return (request_jacobian);
}




template <typename ValType>
bool
MAST::StructuralElementBase::
small_disturbance_surface_pressure_residual(bool request_jacobian,
                                            RealVectorX &f,
                                            RealMatrixX &jac,
                                            const unsigned int side,
                                            MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    libmesh_assert_equal_to(bc.type(), MAST::SMALL_DISTURBANCE_MOTION);
    
    
    MAST::FieldFunction<Real>&
    press_fn = bc.get<MAST::FieldFunction<Real> >("pressure");
    MAST::FieldFunction<ValType>&
    dpress_fn = bc.get<MAST::FieldFunction<ValType> >("dpressure");
    MAST::FieldFunction<typename VectorType<ValType>::return_type>&
    dn_rot_fn = bc.get<MAST::FieldFunction<typename VectorType<ValType>::return_type> >("dnormal");
    
    
    libMesh::FEBase *fe_ptr    = NULL;
    libMesh::QBase  *qrule_ptr = NULL;
    _get_side_fe_and_qrule(get_elem_for_quadrature(),
                           side,
                           &fe_ptr,
                           &qrule_ptr,
                           false);
    std::auto_ptr<libMesh::FEBase> fe(fe_ptr);
    std::auto_ptr<libMesh::QBase>  qrule(qrule_ptr);
    
    
    // Physical location of the quadrature points
    const std::vector<Real> &JxW                    = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint       = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi      = fe->get_phi();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals();
    
    const unsigned int
    n_phi = (unsigned int)phi.size(),
    n1    = 3,
    n2    = 6*n_phi;
    
    RealVectorX phi_vec   = RealVectorX::Zero(n_phi);
    Eigen::Matrix<ValType, Dynamic, 1>
    dn_rot  = Eigen::Matrix<ValType, Dynamic, 1>::Zero(3),
    force   = Eigen::Matrix<ValType, Dynamic, 1>::Zero(2*n1),
    local_f = Eigen::Matrix<ValType, Dynamic, 1>::Zero(n2),
    vec_n2  = Eigen::Matrix<ValType, Dynamic, 1>::Zero(n2);
    
    
    FEMOperatorMatrix Bmat;
    libMesh::Point pt;
    Real press;
    ValType dpress;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location(qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure and deformation information
        press_fn (pt, _time,  press);
        dpress_fn(pt, _time, dpress);
        dn_rot_fn(pt, _time, dn_rot);
        
        //            press = 0.;
        //            dpress = Complex(2./4.*std::real(dn_rot(0)),  2./4./.1*std::imag(utrans(1)));
        //            libMesh::out << q_point[qp](0)
        //            << "  " << std::real(utrans(1))
        //            << "  " << std::imag(utrans(1))
        //            << "  " << std::real(dn_rot(0))
        //            << "  " << std::imag(dn_rot(0))
        //            << "  " << std::real(press)
        //            << "  " << std::imag(press)
        //            << "  " << std::real(dpress)
        //            << "  " << std::imag(dpress) << std::endl;
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) =  ( press * dn_rot(i_dim) + // steady pressure
                             dpress * face_normals[qp](i_dim)); // unsteady pressure
        
        
        Bmat.vector_mult_transpose(vec_n2, force);
        
        local_f -= JxW[qp] * vec_n2;
    }
    
    // now transform to the global system and add
    if (_elem.dim() < 3) {
        transform_vector_to_global_system(local_f, vec_n2);
        MAST::add_to_assembled_vector(f, vec_n2);
    }
    else
        MAST::add_to_assembled_vector(f, local_f);
    
    
    return (request_jacobian);
}







template <typename ValType>
bool
MAST::StructuralElementBase::
small_disturbance_surface_pressure_residual(bool request_jacobian,
                                            RealVectorX &f,
                                            RealMatrixX &jac,
                                            MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements.
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    libmesh_assert_equal_to(bc.type(), MAST::SMALL_DISTURBANCE_MOTION);
    
    MAST::FieldFunction<Real>&
    press_fn = bc.get<MAST::FieldFunction<Real> >("pressure");
    MAST::FieldFunction<ValType>&
    dpress_fn = bc.get<MAST::FieldFunction<ValType> >("dpressure");
    MAST::FieldFunction<typename VectorType<ValType>::return_type>&
    dn_rot_fn = bc.get<MAST::FieldFunction<typename VectorType<ValType>::return_type> >("dnormal");
    
    
    const std::vector<Real> &JxW                 = _fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = _fe->get_phi();
    const unsigned int
    n_phi = (unsigned int)phi.size(),
    n1    = 3,
    n2    = 6*n_phi;
    
    // normal for face integration
    libMesh::Point normal;
    // direction of pressure assumed to be normal (along local z-axis)
    // to the element face for 2D and along local y-axis for 1D element.
    normal(_elem.dim()) = -1.;
    
    RealVectorX phi_vec = RealVectorX::Zero(n_phi);
    Eigen::Matrix<ValType, Dynamic, 1>
    dn_rot  = Eigen::Matrix<ValType, Dynamic, 1>::Zero(3),
    force   = Eigen::Matrix<ValType, Dynamic, 1>::Zero(2*n1),
    local_f = Eigen::Matrix<ValType, Dynamic, 1>::Zero(n2),
    vec_n2  = Eigen::Matrix<ValType, Dynamic, 1>::Zero(n2);
    
    FEMOperatorMatrix Bmat;
    libMesh::Point pt;
    Real press;
    ValType dpress;
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++)
    {
        _local_elem->global_coordinates_location(qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure and deformation information
        press_fn (pt, _time,  press);
        dpress_fn(pt, _time, dpress);
        dn_rot_fn(pt, _time, dn_rot);
        
        //        libMesh::out << std::setw(15) << pt(0)
        //        << std::setw(15) << std::real(press)
        //        << std::setw(15) << std::imag(press)
        //        << std::setw(15) << std::real(dpress)
        //        << std::setw(15) << std::imag(dpress)
        //        << std::setw(15) << std::real(utrans(1))
        //        << std::setw(15) << std::imag(utrans(1))
        //        << std::setw(15) << std::real(dn_rot(0))
        //        << std::setw(15) << std::imag(dn_rot(0)) << std::endl;
        
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = ( press * dn_rot(i_dim) + // steady pressure
                            dpress * normal(i_dim)); // unsteady pressure
        
        Bmat.vector_mult_transpose(vec_n2, force);
        
        local_f -= JxW[qp] * vec_n2;
    }
    
    // now transform to the global system and add
    transform_vector_to_global_system(local_f, vec_n2);
    MAST::add_to_assembled_vector(f, vec_n2);
    
    
    return (request_jacobian);
}









template <typename ValType>
void
MAST::StructuralElementBase::
transform_matrix_to_global_system(const ValType& local_mat,
                                  ValType& global_mat) const {
    
    libmesh_assert_equal_to( local_mat.rows(),  local_mat.cols());
    libmesh_assert_equal_to(global_mat.rows(), global_mat.cols());
    libmesh_assert_equal_to( local_mat.rows(), global_mat.rows());
    
    const unsigned int n_dofs = _fe->n_shape_functions();

    ValType mat(6*n_dofs, 6*n_dofs);
    
    mat.setZero();
    global_mat.setZero();
    
    const RealMatrixX& Tmat = this->local_elem().T_matrix();
    // now initialize the global T matrix
    for (unsigned int i=0; i<n_dofs; i++)
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++) {
                mat(j*n_dofs+i, k*n_dofs+i) = Tmat(j,k); // for u,v,w
                mat((j+3)*n_dofs+i, (k+3)*n_dofs+i) = Tmat(j,k); // for tx,ty,tz
            }
    
    // right multiply with T^tr, and left multiply with T.
    global_mat = mat * local_mat * mat.transpose();
}



template <typename ValType>
void
MAST::StructuralElementBase::
transform_vector_to_local_system(const ValType& global_vec,
                                 ValType& local_vec) const {
    
    libmesh_assert_equal_to( local_vec.size(),  global_vec.size());
    
    const unsigned int n_dofs = _fe->n_shape_functions();
    RealMatrixX mat  = RealMatrixX::Zero(6*n_dofs, 6*n_dofs);
    
    local_vec.setZero();

    const RealMatrixX& Tmat = this->local_elem().T_matrix();
    // now initialize the global T matrix
    for (unsigned int i=0; i<n_dofs; i++)
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++) {
                mat(j*n_dofs+i, k*n_dofs+i) = Tmat(j,k); // for u,v,w
                mat((j+3)*n_dofs+i, (k+3)*n_dofs+i) = Tmat(j,k); // for tx,ty,tz
            }
    
    // left multiply with T^tr
    local_vec = mat.transpose() * global_vec;
}



template <typename ValType>
void
MAST::StructuralElementBase::
transform_vector_to_global_system(const ValType& local_vec,
                                  ValType& global_vec) const {
    
    libmesh_assert_equal_to( local_vec.size(),  global_vec.size());
    
    const unsigned int n_dofs = _fe->n_shape_functions();

    RealMatrixX mat  = RealMatrixX::Zero(6*n_dofs, 6*n_dofs);

    global_vec.setZero();
    
    const RealMatrixX& Tmat = this->local_elem().T_matrix();
    // now initialize the global T matrix
    for (unsigned int i=0; i<n_dofs; i++)
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++) {
                mat(j*n_dofs+i, k*n_dofs+i) = Tmat(j,k); // for u,v,w
                mat((j+3)*n_dofs+i, (k+3)*n_dofs+i) = Tmat(j,k); // for tx,ty,tz
            }
    
    // left multiply with T
    global_vec = mat* local_vec;
}





std::auto_ptr<MAST::StructuralElementBase>
MAST::build_structural_element(MAST::SystemInitialization& sys,
                               const libMesh::Elem& elem,
                               const MAST::ElementPropertyCardBase& p) {
    
    std::auto_ptr<MAST::StructuralElementBase> e;
    
    switch (elem.dim()) {
        case 1:
            e.reset(new MAST::StructuralElement1D(sys, elem, p));
            break;
            
        case 2:
            e.reset(new MAST::StructuralElement2D(sys, elem, p));
            break;
            
        case 3:
            e.reset(new MAST::StructuralElement3D(sys, elem, p));
            break;
            
        default:
            libmesh_error();
            break;
    }
    
    return e;
}



Real
MAST::StructuralElementBase::piston_theory_cp(const unsigned int order,
                                              const Real vel_normal,
                                              const Real a_inf,
                                              const Real gamma,
                                              const Real mach) {
    
    
    Real cp     = 0.0;
    switch (order)
    {
        case 3:
            cp  += (gamma+1.0)/12.0*pow(vel_normal/a_inf,3);
        case 2:
            cp  += (gamma+1.0)/4.0*pow(vel_normal/a_inf,2);
        case 1: {
            cp  += vel_normal/a_inf;
            cp  *= 2.0/pow(mach,2);
        }
            break;
            
        default:
            libmesh_error_msg("Invalid Piston Theory Order: " << order);
            break;
    }
    
    return cp;
}




Real
MAST::StructuralElementBase::piston_theory_dcp_dvn(const unsigned int order,
                                                   const Real vel_normal,
                                                   const Real a_inf,
                                                   const Real gamma,
                                                   const Real mach) {
    
    
    Real dcp_dvn     = 0.0;
    switch (order)
    {
        case 3:
            dcp_dvn  += (gamma+1.0)/4.0*pow(vel_normal,2)/pow(a_inf,3);
        case 2:
            dcp_dvn  += (gamma+1.0)/2.0*vel_normal/pow(a_inf,2);
        case 1: {
            dcp_dvn  += 1./a_inf;
            dcp_dvn  *= 2.0/pow(mach,2);
        }
            break;
            
        default:
            libmesh_error_msg("Invalid Piston Theory Order: " << order);
            break;
    }
    
    return dcp_dvn;
}




// template instantiations
template
void
MAST::StructuralElementBase::transform_matrix_to_global_system<RealMatrixX>
(const RealMatrixX& local_mat,
 RealMatrixX& global_mat) const;


template
void
MAST::StructuralElementBase::transform_vector_to_local_system<RealVectorX>
(const RealVectorX& global_vec,
 RealVectorX& local_vec) const;


template
void
MAST::StructuralElementBase::transform_vector_to_global_system<RealVectorX>
(const RealVectorX& local_vec,
 RealVectorX& global_vec) const;


template
void
MAST::StructuralElementBase::transform_matrix_to_global_system<ComplexMatrixX>
(const ComplexMatrixX& local_mat,
 ComplexMatrixX& global_mat) const;


template
void
MAST::StructuralElementBase::transform_vector_to_local_system<ComplexVectorX>
(const ComplexVectorX& global_vec,
 ComplexVectorX& local_vec) const;


template
void
MAST::StructuralElementBase::transform_vector_to_global_system<ComplexVectorX>
(const ComplexVectorX& local_vec,
 ComplexVectorX& global_vec) const;


template
bool
MAST::StructuralElementBase::side_external_residual<Real>
(bool request_jacobian,
 RealVectorX &f,
 RealMatrixX& jac_xdot,
 RealMatrixX &jac,
 std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> &bc);


template
bool
MAST::StructuralElementBase::side_external_residual<Complex>
(bool request_jacobian,
 RealVectorX &f,
 RealMatrixX& jac_xdot,
 RealMatrixX &jac,
 std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> &bc);


template
bool
MAST::StructuralElementBase::side_external_residual_sensitivity<Real>
(bool request_jacobian,
 RealVectorX &f,
 RealMatrixX& jac_xdot,
 RealMatrixX &jac,
 std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> &bc);


template
bool
MAST::StructuralElementBase::side_external_residual_sensitivity<Complex>
(bool request_jacobian,
 RealVectorX &f,
 RealMatrixX& jac_xdot,
 RealMatrixX &jac,
 std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> &bc);


template
bool
MAST::StructuralElementBase::volume_external_residual<Real>
(bool request_jacobian,
 RealVectorX& f,
 RealMatrixX& jac_xdot,
 RealMatrixX& jac,
 std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);


template
bool
MAST::StructuralElementBase::volume_external_residual<Complex>
(bool request_jacobian,
 RealVectorX& f,
 RealMatrixX& jac_xdot,
 RealMatrixX& jac,
 std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);


template
bool
MAST::StructuralElementBase::volume_external_residual_sensitivity<Real>
(bool request_jacobian,
 RealVectorX& f,
 RealMatrixX& jac_xdot,
 RealMatrixX& jac,
 std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);


template
bool
MAST::StructuralElementBase::volume_external_residual_sensitivity<Complex>
(bool request_jacobian,
 RealVectorX& f,
 RealMatrixX& jac_xdot,
 RealMatrixX& jac,
 std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc);


