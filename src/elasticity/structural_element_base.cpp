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
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_element_1d.h"
#include "elasticity/structural_element_2d.h"
#include "elasticity/solid_element_3d.h"
#include "elasticity/stress_output_base.h"
#include "base/system_initialization.h"
#include "base/boundary_condition_base.h"
#include "base/nonlinear_system.h"
#include "base/assembly_base.h"
#include "base/field_function_base.h"
#include "property_cards/element_property_card_base.h"
#include "mesh/geom_elem.h"
#include "mesh/fe_base.h"
#include "numerics/fem_operator_matrix.h"
#include "numerics/utility.h"


MAST::StructuralElementBase::StructuralElementBase(MAST::SystemInitialization& sys,
                                                   const MAST::GeomElem& elem,
                                                   const MAST::ElementPropertyCardBase& p):
MAST::ElementBase(sys, elem),
follower_forces   (false),
_property         (p),
_incompatible_sol (nullptr) {
    
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
MAST::StructuralElementBase::set_perturbed_solution(const RealVectorX& vec,
                                                    bool if_sens) {
    
    // convert the vector to the local element coordinate system
    if (!if_sens) {
        if (_elem.dim() == 3)
            _local_delta_sol = vec;
        else {
            _local_delta_sol = RealVectorX::Zero(vec.size());
            this->transform_vector_to_local_system(vec, _local_delta_sol);
        }
    }
    else {
        
        // set the element solution sensitivity.
        if (_elem.dim() == 3)
            _local_delta_sol_sens = vec;
        else {
            _local_delta_sol_sens = RealVectorX::Zero(vec.size());
            this->transform_vector_to_local_system(vec, _local_delta_sol_sens);
        }
    }
    
    MAST::ElementBase::set_perturbed_solution(vec, if_sens);
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
MAST::StructuralElementBase::set_perturbed_velocity(const RealVectorX& vec,
                                                    bool if_sens) {
    
    if (!if_sens) {
        if (_elem.dim() == 3)
            _local_delta_vel = vec;
        else {
            _local_delta_vel = RealVectorX::Zero(vec.size());
            this->transform_vector_to_local_system(vec, _local_delta_vel);
        }
    }
    else {
        
        if (_elem.dim() == 3)
            _local_delta_vel_sens = vec;
        else {
            _local_delta_vel_sens = RealVectorX::Zero(vec.size());
            this->transform_vector_to_local_system(vec, _local_delta_vel_sens);
        }
    }
    
    MAST::ElementBase::set_perturbed_velocity(vec, if_sens);
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



void
MAST::StructuralElementBase::set_perturbed_acceleration(const RealVectorX& vec,
                                                        bool if_sens) {
    
    if (!if_sens) {
        if (_elem.dim() == 3)
            _local_delta_accel = vec;
        else {
            _local_delta_accel = RealVectorX::Zero(vec.size());
            this->transform_vector_to_local_system(vec, _local_delta_accel);
        }
    }
    else {
        
        if (_elem.dim() == 3)
            _local_delta_accel_sens = vec;
        else {
            _local_delta_accel_sens = RealVectorX::Zero(vec.size());
            this->transform_vector_to_local_system(vec, _local_delta_accel_sens);
        }
    }
    
    MAST::ElementBase::set_perturbed_acceleration(vec, if_sens);
}



bool
MAST::StructuralElementBase::linearized_internal_residual (bool request_jacobian,
                                                           RealVectorX& f,
                                                           RealMatrixX& jac) {
    
    
    // It is assumed that the structural elements implement the Jacobians.
    // Hence, the residual for the linearized problem will be calculated
    // using the Jacobian.
    
    const unsigned int
    n_phi    = (unsigned int)f.size(),
    n2       =6*n_phi;

    RealVectorX
    v = RealVectorX::Zero(n2);
    
    RealMatrixX
    m = RealMatrixX::Zero(n2, n2);
    
    
    this->internal_residual(true, v, m);
    
    f += m * _local_delta_sol;
    
    if (request_jacobian)
        jac += m;
    
    return request_jacobian;
}



bool
MAST::StructuralElementBase::inertial_residual (bool request_jacobian,
                                                RealVectorX& f,
                                                RealMatrixX& jac_xddot,
                                                RealMatrixX& jac_xdot,
                                                RealMatrixX& jac) {
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real>& JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz     = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    
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
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
    mat_inertia  = _property.inertia_matrix(*this);
    
    MAST::FEMOperatorMatrix Bmat;
    
    if (_property.if_diagonal_mass_matrix()) {
        
        (*mat_inertia)(xyz[0], _time, material_mat);
        
        Real vol = 0.;
        const unsigned int nshp = fe->n_shape_functions();
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
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            
            (*mat_inertia)(xyz[qp], _time, material_mat);
            
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
MAST::StructuralElementBase::
linearized_inertial_residual (bool request_jacobian,
                              RealVectorX& f,
                              RealMatrixX& jac_xddot,
                              RealMatrixX& jac_xdot,
                              RealMatrixX& jac) {
    
    // It is assumed that the structural elements implement the Jacobians.
    // Hence, the residual for the linearized problem will be calculated
    // using the Jacobian.
    
    const unsigned int
    n_phi    = (unsigned int)f.size(),
    n2       =6*n_phi;
    
    RealVectorX
    v = RealVectorX::Zero(n2);
    
    RealMatrixX
    m_xddot = RealMatrixX::Zero(n2, n2),
    m_xdot  = RealMatrixX::Zero(n2, n2),
    m_x     = RealMatrixX::Zero(n2, n2);
    
    
    this->inertial_residual(true, v, m_xddot, m_xdot, m_x);
    
    f += (m_xddot * _local_delta_accel +
          m_xdot  * _local_delta_vel   +
          m_x     * _local_delta_sol);
    
    if (request_jacobian) {
        
        jac_xddot += m_xddot;
        jac_xdot  +=  m_xdot;
        jac       +=     m_x;
    }
    
    return request_jacobian;
}




bool
MAST::StructuralElementBase::inertial_residual_sensitivity (const MAST::FunctionBase& p,
                                                            bool request_jacobian,
                                                            RealVectorX& f,
                                                            RealMatrixX& jac_xddot,
                                                            RealMatrixX& jac_xdot,
                                                            RealMatrixX& jac) {
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real>& JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz     = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    
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
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
    mat_inertia  = _property.inertia_matrix(*this);
    
    MAST::FEMOperatorMatrix Bmat;
    
    if (_property.if_diagonal_mass_matrix()) {
        
        // as an approximation, get matrix at the first quadrature point
        mat_inertia->derivative(p, xyz[0], _time, material_mat);
        
        Real vol = 0.;
        const unsigned int nshp = fe->n_shape_functions();
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
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            
            mat_inertia->derivative(p, xyz[qp], _time, material_mat);
            
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




void
MAST::StructuralElementBase::
inertial_residual_boundary_velocity (const MAST::FunctionBase& p,
                                     const unsigned int s,
                                     const MAST::FieldFunction<RealVectorX>& vel_f,
                                     bool request_jacobian,
                                     RealVectorX& f,
                                     RealMatrixX& jac_xddot,
                                     RealMatrixX& jac_xdot,
                                     RealMatrixX& jac) {
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    std::vector<Real> JxW_Vn                        = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz          = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi      = fe->get_phi();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals_for_local_coordinate();

    const unsigned int
    n_phi    = (unsigned int)phi.size(),
    n_vars   = _system.system().n_vars(),
    n1       =6,
    n2       =6*n_phi,
    dim      =_elem.dim();
    
    RealMatrixX
    material_mat,
    mat1_n1n2     = RealMatrixX::Zero(n1, n2),
    mat2_n2n2     = RealMatrixX::Zero(n2, n2),
    local_jac     = RealMatrixX::Zero(n2, n2);
    RealVectorX
    phi_vec    = RealVectorX::Zero(n_phi),
    vec1_n1    = RealVectorX::Zero(n1),
    vec2_n2    = RealVectorX::Zero(n2),
    local_f    = RealVectorX::Zero(n2),
    vel        = RealVectorX::Zero(dim);
    
    Real
    vn  = 0.;
    
    // modify the JxW_Vn by multiplying the normal velocity to it
    for (unsigned int qp=0; qp<JxW_Vn.size(); qp++) {
        
        vel_f.derivative(p, xyz[qp], _time, vel);
        vn = 0.;
        for (unsigned int i=0; i<dim; i++)
            vn += vel(i)*face_normals[qp](i);
        JxW_Vn[qp] *= vn;
    }
    

    std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
    mat_inertia  = _property.inertia_matrix(*this);
    
    MAST::FEMOperatorMatrix Bmat;
    
    if (_property.if_diagonal_mass_matrix()) {
        
        (*mat_inertia)(xyz[0], _time, material_mat);
        
        Real vol = 0.;
        const unsigned int nshp = fe->n_shape_functions();
        for (unsigned int i=0; i<JxW_Vn.size(); i++)
            vol += JxW_Vn[i];
        vol /= (1.* nshp);
        for (unsigned int i_var=0; i_var<6; i_var++)
            for (unsigned int i=0; i<nshp; i++)
                local_jac(i_var*nshp+i, i_var*nshp+i) =
                vol*material_mat(i_var, i_var);
        
        local_f =  local_jac * _local_accel;
    }
    else {
        
        for (unsigned int qp=0; qp<JxW_Vn.size(); qp++) {
            
            (*mat_inertia)(xyz[qp], _time, material_mat);
            
            // now set the shape function values
            for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
                phi_vec(i_nd) = phi[i_nd][qp];
            
            Bmat.reinit(n_vars, phi_vec);
            
            Bmat.left_multiply(mat1_n1n2, material_mat);
            
            vec1_n1 = mat1_n1n2 * _local_accel;
            Bmat.vector_mult_transpose(vec2_n2, vec1_n1);
            
            local_f += JxW_Vn[qp] * vec2_n2;
            
            if (request_jacobian) {
                
                Bmat.right_multiply_transpose(mat2_n2n2,
                                              mat1_n1n2);
                local_jac += JxW_Vn[qp]*mat2_n2n2;
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
}




bool
MAST::StructuralElementBase::
side_external_residual(bool request_jacobian,
                       RealVectorX& f,
                       RealMatrixX& jac_xdot,
                       RealMatrixX& jac,
                       std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>> loads;
    _elem.external_side_loads_for_quadrature_elem(bc, loads);
    
    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>>::const_iterator
    it   = loads.begin(),
    end  = loads.end();
    
    for ( ; it != end; it++) {
        
        std::vector<MAST::BoundaryConditionBase*>::const_iterator
        bc_it  = it->second.begin(),
        bc_end = it->second.end();
        
        for ( ; bc_it != bc_end; bc_it++) {
            
            // apply all the types of loading
            switch ((*bc_it)->type()) {
                    
                case MAST::SURFACE_PRESSURE:
                    surface_pressure_residual(request_jacobian,
                                              f, jac,
                                              it->first,
                                              **bc_it);
                    break;


                case MAST::SURFACE_TRACTION:
                    surface_traction_residual(request_jacobian,
                                              f, jac,
                                              it->first,
                                              **bc_it);
                    break;

                    
                case MAST::SURFACE_TRACTION_SHIFTED_BOUNDARY:
                    surface_traction_residual_shifted_boundary(request_jacobian,
                                                               f, jac,
                                                               it->first,
                                                               **bc_it);
                    break;
                    
                    
                case MAST::PISTON_THEORY:
                    piston_theory_residual(request_jacobian,
                                           f,
                                           jac_xdot,
                                           jac,
                                           it->first,
                                           **bc_it);
                    break;
                    
                    
                case MAST::BOUNDARY_VELOCITY:
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
    return request_jacobian;
}

bool
MAST::StructuralElementBase::
linearized_side_external_residual
(bool request_jacobian,
 RealVectorX& f,
 RealMatrixX& jac_xdot,
 RealMatrixX& jac,
 std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {

    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>> loads;
    _elem.external_side_loads_for_quadrature_elem(bc, loads);
    
    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>>::const_iterator
    it   = loads.begin(),
    end  = loads.end();
    
    for ( ; it != end; it++) {
        
        std::vector<MAST::BoundaryConditionBase*>::const_iterator
        bc_it  = it->second.begin(),
        bc_end = it->second.end();
        
        for ( ; bc_it != bc_end; bc_it++) {
            
            // apply all the types of loading
            switch ((*bc_it)->type()) {

                case MAST::BOUNDARY_VELOCITY:
                case MAST::DIRICHLET:
                    // nothing to be done here
                    break;
                    
                case MAST::SURFACE_PRESSURE:
                case MAST::PISTON_THEORY:
                default:
                    // not implemented yet
                    libmesh_error();
                    break;
            }
        }
    }
    return request_jacobian;
}



bool
MAST::StructuralElementBase::
linearized_frequency_domain_side_external_residual
(bool request_jacobian,
 ComplexVectorX& f,
 ComplexMatrixX& jac,
 std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>> loads;
    _elem.external_side_loads_for_quadrature_elem(bc, loads);
    
    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>>::const_iterator
    it   = loads.begin(),
    end  = loads.end();
    
    for ( ; it != end; it++) {
        
        std::vector<MAST::BoundaryConditionBase*>::const_iterator
        bc_it  = it->second.begin(),
        bc_end = it->second.end();
        
        for ( ; bc_it != bc_end; bc_it++) {
            
            // apply all the types of loading
            switch ((*bc_it)->type()) {
                    
                case MAST::SURFACE_PRESSURE:
                    
                    linearized_frequency_domain_surface_pressure_residual
                    (request_jacobian,
                     f, jac,
                     it->first,
                     **bc_it);
                    break;
                    
                    
                case MAST::BOUNDARY_VELOCITY:
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
    return request_jacobian;
}




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
    
    libMesh::subdomain_id_type sid = _elem.get_reference_elem().subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::SURFACE_PRESSURE:
                surface_pressure_residual(request_jacobian,
                                          f, jac,
                                          *it.first->second);
                break;

            case MAST::PISTON_THEORY:
                piston_theory_residual(request_jacobian,
                                       f,
                                       jac_xdot,
                                       jac,
                                       *it.first->second);
                break;

            case MAST::TEMPERATURE:
                thermal_residual(request_jacobian,
                                 f, jac,
                                 *it.first->second);
                break;
                

            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    
    return request_jacobian;
}




bool
MAST::StructuralElementBase::
linearized_volume_external_residual (bool request_jacobian,
                                     RealVectorX& f,
                                     RealMatrixX& jac_xdot,
                                     RealMatrixX& jac,
                                     std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    // iterate over the boundary ids given in the provided force map
    std::pair<std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator,
    std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    
    libMesh::subdomain_id_type sid = _elem.get_reference_elem().subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::SURFACE_PRESSURE:
                linearized_surface_pressure_residual(request_jacobian,
                                                     f, jac,
                                                     *it.first->second);
                break;
                
            case MAST::TEMPERATURE: {
                
                const unsigned int
                n = (unsigned int) f.size();
                
                RealVectorX
                local_f  =  RealVectorX::Zero(n);
                RealMatrixX
                mat      =  RealMatrixX::Zero(n, n);
                
                // this accounts for the perturbation of displacement.
                // Perturbation in temperature is not currently handled
                thermal_residual(true, local_f, mat, *it.first->second);
                f  +=  mat*_delta_sol;
            }
                break;
                
            case MAST::PISTON_THEORY:
            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    
    return request_jacobian;
}




bool
MAST::StructuralElementBase::
linearized_frequency_domain_volume_external_residual
(bool request_jacobian,
 ComplexVectorX& f,
 ComplexMatrixX& jac,
 std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    // iterate over the boundary ids given in the provided force map
    std::pair<std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator,
    std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    
    libMesh::subdomain_id_type sid = _elem.get_reference_elem().subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::SURFACE_PRESSURE:
                
                linearized_frequency_domain_surface_pressure_residual
                (request_jacobian,
                 f, jac,
                 *it.first->second);
                break;
                
            case MAST::TEMPERATURE:
                
            case MAST::PISTON_THEORY:
            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    
    return request_jacobian;
}







bool
MAST::StructuralElementBase::
side_external_residual_sensitivity(const MAST::FunctionBase& p,
                                   bool request_jacobian,
                                   RealVectorX& f,
                                   RealMatrixX& jac_xdot,
                                   RealMatrixX& jac,
                                   std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>> loads;
    _elem.external_side_loads_for_quadrature_elem(bc, loads);
    
    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>>::const_iterator
    it   = loads.begin(),
    end  = loads.end();
    
    for ( ; it != end; it++) {
        
        std::vector<MAST::BoundaryConditionBase*>::const_iterator
        bc_it  = it->second.begin(),
        bc_end = it->second.end();
        
        for ( ; bc_it != bc_end; bc_it++) {
            
            // apply all the types of loading
            switch ((*bc_it)->type()) {

                case MAST::SURFACE_PRESSURE:
                    surface_pressure_residual_sensitivity(p,
                                                          request_jacobian,
                                                          f, jac,
                                                          it->first,
                                                          **bc_it);
                    break;
                    
                    
                case MAST::SURFACE_TRACTION:
                    surface_traction_residual_sensitivity(p,
                                                          request_jacobian,
                                                          f, jac,
                                                          it->first,
                                                          **bc_it);
                    break;
                    
                    
                case MAST::SURFACE_TRACTION_SHIFTED_BOUNDARY:
                    surface_traction_residual_shifted_boundary_sensitivity(p,
                                                                           request_jacobian,
                                                                           f, jac,
                                                                           it->first,
                                                                           **bc_it);
                    break;
                    
                    
                case MAST::PISTON_THEORY:
                    piston_theory_residual_sensitivity(p,
                                                       request_jacobian,
                                                       f,
                                                       jac_xdot,
                                                       jac,
                                                       it->first,
                                                       **bc_it);
                    break;
                    
                    
                case MAST::BOUNDARY_VELOCITY:
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
    return request_jacobian;
}



bool
MAST::StructuralElementBase::
volume_external_residual_sensitivity (const MAST::FunctionBase& p,
                                      bool request_jacobian,
                                      RealVectorX& f,
                                      RealMatrixX& jac_xdot,
                                      RealMatrixX& jac,
                                      std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    // iterate over the boundary ids given in the provided force map
    std::pair<std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator,
    std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    
    libMesh::subdomain_id_type sid = _elem.get_reference_elem().subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::SURFACE_PRESSURE:
                surface_pressure_residual_sensitivity(p,
                                                      request_jacobian,
                                                      f, jac,
                                                      *it.first->second);
                break;
                
            case MAST::PISTON_THEORY:
                piston_theory_residual_sensitivity(p,
                                                   request_jacobian,
                                                   f,
                                                   jac_xdot,
                                                   jac,
                                                   *it.first->second);
                break;
                
            case MAST::TEMPERATURE:
                thermal_residual_sensitivity(p,
                                             request_jacobian,
                                             f, jac,
                                             *it.first->second);
                break;
                
            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    
    return request_jacobian;
}



void
MAST::StructuralElementBase::
volume_external_residual_boundary_velocity(const MAST::FunctionBase& p,
                                           const unsigned int s,
                                           const MAST::FieldFunction<RealVectorX>& vel_f,
                                           std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc,
                                           bool request_jacobian,
                                           RealVectorX& f,
                                           RealMatrixX& jac) {
    
    // iterate over the boundary ids given in the provided force map
    std::pair<std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator,
    std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    
    libMesh::subdomain_id_type sid = _elem.get_reference_elem().subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::SURFACE_PRESSURE:
                surface_pressure_boundary_velocity(p,
                                                   s,
                                                   vel_f,
                                                   *it.first->second,
                                                   request_jacobian,
                                                   f,
                                                   jac);
                break;
                
            case MAST::TEMPERATURE:
                thermal_residual_boundary_velocity(p,
                                                   s,
                                                   vel_f,
                                                   *it.first->second,
                                                   request_jacobian,
                                                   f,
                                                   jac);
                break;
                
            case MAST::PISTON_THEORY:
            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
}



bool
MAST::StructuralElementBase::
surface_pressure_residual(bool request_jacobian,
                          RealVectorX &f,
                          RealMatrixX &jac,
                          MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));
    
    const std::vector<Real> &JxW                = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint   = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi  = fe->get_phi();
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
    
    RealVectorX
    phi_vec  = RealVectorX::Zero(n_phi),
    force    = RealVectorX::Zero(2*n1),
    local_f  = RealVectorX::Zero(n2),
    vec_n2   = RealVectorX::Zero(n2);
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure value
        func(qpoint[qp], _time, press);
        
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
surface_pressure_residual_sensitivity(const MAST::FunctionBase& p,
                                      bool request_jacobian,
                                      RealVectorX &f,
                                      RealMatrixX &jac,
                                      MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real> &JxW                = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint   = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi  = fe->get_phi();
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
    
    RealVectorX
    phi_vec  = RealVectorX::Zero(n_phi),
    force    = RealVectorX::Zero(2*n1),
    local_f  = RealVectorX::Zero(n2),
    vec_n2   = RealVectorX::Zero(n2);
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure value
        func.derivative(p, qpoint[qp], _time, press);
        
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




void
MAST::StructuralElementBase::
surface_pressure_boundary_velocity(const MAST::FunctionBase& p,
                                   const unsigned int s,
                                   const MAST::FieldFunction<RealVectorX>& vel_f,
                                   MAST::BoundaryConditionBase& bc,
                                   bool request_jacobian,
                                   RealVectorX& f,
                                   RealMatrixX& jac) {
    
    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    std::vector<Real> JxW_Vn                        = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz          = fe->get_xyz();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals_for_local_coordinate();
    const std::vector<std::vector<Real> >& phi      = fe->get_phi();

    const unsigned int
    n_phi = (unsigned int)phi.size(),
    n1    = 3,
    n2    = 6*n_phi,
    dim   = _elem.dim();
    
    
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
    
    RealVectorX
    phi_vec  = RealVectorX::Zero(n_phi),
    force    = RealVectorX::Zero(2*n1),
    local_f  = RealVectorX::Zero(n2),
    vec_n2   = RealVectorX::Zero(n2),
    vel      = RealVectorX::Zero(dim);

    Real
    vn  = 0.;
    
    
    // modify the JxW_Vn by multiplying the normal velocity to it
    for (unsigned int qp=0; qp<JxW_Vn.size(); qp++) {
        
        vel_f.derivative(p, xyz[qp], _time, vel);
        vn = 0.;
        for (unsigned int i=0; i<dim; i++)
            vn += vel(i)*face_normals[qp](i);
        JxW_Vn[qp] *= vn;
    }

    
    for (unsigned int qp=0; qp<xyz.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure value
        func(xyz[qp], _time, press);
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = press * normal(i_dim);
        
        Bmat.vector_mult_transpose(vec_n2, force);
        
        local_f += JxW_Vn[qp] * vec_n2;
    }
    
    // now transform to the global system and add
    transform_vector_to_global_system(local_f, vec_n2);
    f -= vec_n2;
}





bool
MAST::StructuralElementBase::
linearized_surface_pressure_residual(bool request_jacobian,
                                     RealVectorX &f,
                                     RealMatrixX &jac,
                                     MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements.
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    libmesh_assert_equal_to(bc.type(), MAST::SURFACE_PRESSURE);
    
    MAST::FieldFunction<Real>&
    press_fn  = bc.get<MAST::FieldFunction<Real> >("pressure");
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = fe->get_phi();
    const unsigned int
    n_phi = (unsigned int)phi.size(),
    n1    = 3,
    n2    = 6*n_phi;
    
    // normal for face integration
    libMesh::Point normal;
    // direction of pressure assumed to be normal (along local z-axis)
    // to the element face for 2D and along local y-axis for 1D element.
    normal(_elem.dim()) = -1.;
    
    RealVectorX
    phi_vec = RealVectorX::Zero(n_phi),
    w       = RealVectorX::Zero(3),
    dn_rot  = RealVectorX::Zero(3),
    force   = RealVectorX::Zero(2*n1),
    local_f = RealVectorX::Zero(n2),
    vec_n2  = RealVectorX::Zero(n2);
    
    FEMOperatorMatrix Bmat;
    Real
    press   = 0.,
    dpress  = 0.;
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure and deformation information
        press_fn             (qpoint[qp], _time,  press);
        press_fn.perturbation(qpoint[qp], _time, dpress);
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = ( press * dn_rot(i_dim) + // steady pressure
                            dpress * normal(i_dim)); // unsteady pressure
        
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
linearized_frequency_domain_surface_pressure_residual
(bool request_jacobian,
 ComplexVectorX &f,
 ComplexMatrixX &jac,
 MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements.
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    libmesh_assert_equal_to(bc.type(), MAST::SURFACE_PRESSURE);
    
    MAST::FieldFunction<Real>&
    press_fn   = bc.get<MAST::FieldFunction<Real> >("pressure");
    MAST::FieldFunction<Complex>&
    dpress_fn  = bc.get<MAST::FieldFunction<Complex> >("frequency_domain_pressure");
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = fe->get_phi();
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
    ComplexVectorX
    w       = ComplexVectorX::Zero(3),
    dn_rot  = ComplexVectorX::Zero(3),
    force   = ComplexVectorX::Zero(2*n1),
    local_f = ComplexVectorX::Zero(n2),
    vec_n2  = ComplexVectorX::Zero(n2);
    
    FEMOperatorMatrix Bmat;
    Real              press  = 0.;
    Complex           dpress = 0.;
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure and deformation information
        press_fn (qpoint[qp], _time,  press);
        dpress_fn(qpoint[qp], _time, dpress);
        //dn_rot_fn.freq_domain_motion(pt, normal, w, dn_rot);
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = ( press * dn_rot(i_dim) + // steady pressure
                            dpress * normal(i_dim)); // unsteady pressure
        
        Bmat.vector_mult_transpose(vec_n2, force);
        
        local_f += JxW[qp] * vec_n2;
    }
    
    // now transform to the global system and add
    transform_vector_to_global_system(local_f, vec_n2);
    f -= vec_n2;
    
    
    return (request_jacobian);
}



template <typename ValType>
void
MAST::StructuralElementBase::
transform_matrix_to_global_system(const ValType& local_mat,
                                  ValType& global_mat) const {
    
    // make sure this is only called for 1D and 2D elements
    libmesh_assert_less(_elem.dim(), 3);
    libmesh_assert_equal_to( local_mat.rows(),  local_mat.cols());
    libmesh_assert_equal_to(global_mat.rows(), global_mat.cols());
    libmesh_assert_equal_to( local_mat.rows(), global_mat.rows());
    
    if (_elem.use_local_elem()) {
        
        // assuming same shape function for all 6 dofs.
        const unsigned int n_dofs = local_mat.rows()/6;
        
        ValType mat(6*n_dofs, 6*n_dofs);
        
        mat.setZero();
        global_mat.setZero();
        
        const RealMatrixX& Tmat = _elem.T_matrix();
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
    else
        global_mat = local_mat;
}



template <typename ValType>
void
MAST::StructuralElementBase::
transform_vector_to_local_system(const ValType& global_vec,
                                 ValType& local_vec) const {
    
    // make sure this is only called for 1D and 2D elements
    libmesh_assert_less(_elem.dim(), 3);
    libmesh_assert_equal_to( local_vec.size(),  global_vec.size());
    
    if (_elem.use_local_elem()) {

        // assuming same shape function for all 6 dofs.
        const unsigned int n_dofs = global_vec.size()/6;
        RealMatrixX mat  = RealMatrixX::Zero(6*n_dofs, 6*n_dofs);
        
        local_vec.setZero();
        
        const RealMatrixX& Tmat = _elem.T_matrix();
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
    else
        local_vec = global_vec;
}



template <typename ValType>
void
MAST::StructuralElementBase::
transform_vector_to_global_system(const ValType& local_vec,
                                  ValType& global_vec) const {
    
    // make sure this is only called for 1D and 2D elements
    libmesh_assert_less(_elem.dim(), 3);
    libmesh_assert_equal_to( local_vec.size(),  global_vec.size());
    
    if (_elem.use_local_elem()) {
        
        // assuming same shape function for all 6 dofs.
        const unsigned int n_dofs = local_vec.size()/6;
        
        RealMatrixX mat  = RealMatrixX::Zero(6*n_dofs, 6*n_dofs);
        
        global_vec.setZero();
        
        const RealMatrixX& Tmat = _elem.T_matrix();
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
    else
        global_vec = local_vec;
}





std::unique_ptr<MAST::StructuralElementBase>
MAST::build_structural_element(MAST::SystemInitialization& sys,
                               const MAST::GeomElem& elem,
                               const MAST::ElementPropertyCardBase& p) {
    
    std::unique_ptr<MAST::StructuralElementBase> e;
    
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
                                              const Real vel_U,
                                              const Real gamma,
                                              const Real mach) {
    
    
    Real cp     = 0.0;
    switch (order)
    {
        case 3:
            cp  += (gamma+1.0)/12.0*mach*mach*pow(vel_U,3);
        case 2:
            cp  += (gamma+1.0)/4.0*mach*pow(vel_U,2);
        case 1: {
            cp  += vel_U;
            cp  *= 2.0/pow(mach*mach-1.,0.5);
        }
            break;
            
        default:
            libmesh_error_msg("Invalid Piston Theory Order: " << order);
            break;
    }
    
    return cp;
}




Real
MAST::StructuralElementBase::piston_theory_dcp_dv(const unsigned int order,
                                                  const Real vel_U,
                                                  const Real gamma,
                                                  const Real mach) {
    
    
    Real dcp_dvn     = 0.0;
    switch (order)
    {
        case 3:
            dcp_dvn  += (gamma+1.0)/4.0*mach*mach*pow(vel_U,2);
        case 2:
            dcp_dvn  += (gamma+1.0)/2.0*mach*vel_U;
        case 1: {
            dcp_dvn  += 1.;
            dcp_dvn  *= 2.0/pow(mach*mach-1.,.5);
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



