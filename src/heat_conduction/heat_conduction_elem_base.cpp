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
#include "heat_conduction/heat_conduction_elem_base.h"
#include "numerics/fem_operator_matrix.h"
#include "base/system_initialization.h"
#include "base/field_function_base.h"
#include "base/parameter.h"
#include "base/boundary_condition_base.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"
#include "base/assembly_base.h"
#include "mesh/fe_base.h"
#include "mesh/geom_elem.h"
#include "property_cards/element_property_card_base.h"


MAST::HeatConductionElementBase::
HeatConductionElementBase(MAST::SystemInitialization&          sys,
                          const MAST::GeomElem&                 elem,
                          const MAST::ElementPropertyCardBase& p):
MAST::ElementBase (sys, elem),
_property         (p) {

}



MAST::HeatConductionElementBase::~HeatConductionElementBase() {
    
}




void
MAST::HeatConductionElementBase::internal_residual (bool request_jacobian,
                                                    RealVectorX& f,
                                                    RealMatrixX& jac) {
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(true, false));

    const std::vector<Real>& JxW           = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = fe->get_xyz();
    const unsigned int
    n_phi  = fe->n_shape_functions(),
    dim    = _elem.dim();
    
    RealMatrixX
    material_mat   = RealMatrixX::Zero(dim, dim),
    dmaterial_mat  = RealMatrixX::Zero(dim, dim), // for calculation of Jac when k is temp. dep.
    mat_n2n2       = RealMatrixX::Zero(n_phi, n_phi);
    RealVectorX
    vec1     = RealVectorX::Zero(1),
    vec2_n2  = RealVectorX::Zero(n_phi),
    flux     = RealVectorX::Zero(dim);
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX> > conductance =
    _property.thermal_conductance_matrix(*this);
    
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix Bmat; // for calculation of Jac when k is temp. dep.

    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {

        // initialize the Bmat operator for this term
        _initialize_mass_fem_operator(qp, *fe, Bmat);
        Bmat.right_multiply(vec1, _sol);

        if (_active_sol_function)
            dynamic_cast<MAST::MeshFieldFunction*>
            (_active_sol_function)->set_element_quadrature_point_solution(vec1);
        
        (*conductance)(xyz[qp], _time, material_mat);

        _initialize_fem_gradient_operator(qp, dim, *fe, dBmat);
        
        // calculate the flux for each dimension and add its weighted
        // component to the residual
        flux.setZero();
        for (unsigned int j=0; j<dim; j++) {
            dBmat[j].right_multiply(vec1, _sol);        // dT_dxj
            
            for (unsigned int i=0; i<dim; i++)
                flux(i) += vec1(0) * material_mat(i,j); // q_i = k_ij dT_dxj
        }

        // now add to the residual vector
        for (unsigned int i=0; i<dim; i++) {
            vec1(0)  = flux(i);
            dBmat[i].vector_mult_transpose(vec2_n2, vec1);
            f += JxW[qp] * vec2_n2;
        }

        
        if (request_jacobian) {
            
            // Jacobian contribution from int_omega dB_dxi^T k_ij dB_dxj
            for (unsigned int i=0; i<dim; i++)
                for (unsigned int j=0; j<dim; j++) {
                    
                    dBmat[i].right_multiply_transpose(mat_n2n2, dBmat[j]);
                    jac += JxW[qp] * material_mat(i,j) * mat_n2n2;
                }
            
            // Jacobian contribution from int_omega dB_dxi dT_dxj dk_ij/dT B
            if (_active_sol_function) {
                // get derivative of the conductance matrix wrt temperature
                conductance->derivative(*_active_sol_function,
                                        xyz[qp],
                                        _time, dmaterial_mat);
                
                for (unsigned int j=0; j<dim; j++) {
                    dBmat[j].right_multiply(vec1, _sol);  // dT_dxj

                    for (unsigned int i=0; i<dim; i++)
                        if (dmaterial_mat(i,j) != 0.) { // no need to process for zero terms
                            // dB_dxi^T B
                            dBmat[i].right_multiply_transpose(mat_n2n2, Bmat);
                            // dB_dxi^T (dT_dxj dk_ij/dT) B
                            jac += JxW[qp] * vec1(0) * dmaterial_mat(i,j) * mat_n2n2;
                        }
                }
            }
        }
    }
    
    if (_active_sol_function)
        dynamic_cast<MAST::MeshFieldFunction*>
        (_active_sol_function)->clear_element_quadrature_point_solution();
}





void
MAST::HeatConductionElementBase::velocity_residual (bool request_jacobian,
                                                    RealVectorX& f,
                                                    RealMatrixX& jac_xdot,
                                                    RealMatrixX& jac) {
    MAST::FEMOperatorMatrix Bmat;
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real>& JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz       = fe->get_xyz();
    
    const unsigned int
    n_phi      = fe->n_shape_functions(),
    dim        = _elem.dim();
    
    RealMatrixX
    material_mat    = RealMatrixX::Zero(dim, dim),
    mat_n2n2        = RealMatrixX::Zero(n_phi, n_phi);
    RealVectorX
    vec1    = RealVectorX::Zero(1),
    vec2_n2 = RealVectorX::Zero(n_phi),
    local_f = RealVectorX::Zero(n_phi);
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX> > capacitance =
    _property.thermal_capacitance_matrix(*this);
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        _initialize_mass_fem_operator(qp, *fe, Bmat);
        Bmat.right_multiply(vec1, _sol);               //  B * T
        
        if (_active_sol_function)
            dynamic_cast<MAST::MeshFieldFunction*>
            (_active_sol_function)->set_element_quadrature_point_solution(vec1);
        
        (*capacitance)(xyz[qp], _time, material_mat);
        
        Bmat.right_multiply(vec1, _vel);               //  B * T_dot
        Bmat.vector_mult_transpose(vec2_n2, vec1);     //  B^T * B * T_dot
        
        local_f    += JxW[qp] * material_mat(0,0) * vec2_n2; // (rho*cp)*JxW B^T B T_dot
        
        if (request_jacobian || _property.if_diagonal_mass_matrix()) {
            
            Bmat.right_multiply_transpose(mat_n2n2, Bmat);  // B^T B
            jac_xdot += JxW[qp] * material_mat(0,0) * mat_n2n2;  // B^T B * JxW (rho*cp)
            
            // Jacobian contribution from int_omega B T d(rho*cp)/dT B
            if (_active_sol_function) {
                // get derivative of the conductance matrix wrt temperature
                capacitance->derivative(*_active_sol_function,
                                        xyz[qp],
                                        _time, material_mat);
                
                if (material_mat(0,0) != 0.) { // no need to process for zero terms
                    
                    // B^T (T d(rho cp)/dT) B
                    jac += JxW[qp] * vec1(0) * material_mat(0,0) * mat_n2n2;
                }
            }
        }
    }
    
    // diagonalize the matrix and compute the residual based on that.
    if (_property.if_diagonal_mass_matrix()) {

        for (unsigned int i=0; i<jac_xdot.rows(); i++) {
            
            Real a = jac_xdot.row(i).sum();
            jac_xdot.row(i).setZero();
            jac_xdot(i,i) = a;
        }

        f +=  jac_xdot * _vel;

        // if the jacobian was not requested, then zero the matrix.
        if (!request_jacobian) {
            jac_xdot.setZero();
            jac.setZero();
        }
    }
    else
        f += local_f;

    if (_active_sol_function)
        dynamic_cast<MAST::MeshFieldFunction*>
        (_active_sol_function)->clear_element_quadrature_point_solution();
}




void
MAST::HeatConductionElementBase::
side_external_residual (bool request_jacobian,
                        RealVectorX& f,
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
                case MAST::HEAT_FLUX:
                    surface_flux_residual(request_jacobian,
                                          f, jac,
                                          it->first,
                                          **bc_it);
                    break;
                    
                case MAST::CONVECTION_HEAT_FLUX:
                    surface_convection_residual(request_jacobian,
                                                f, jac,
                                                it->first,
                                                **bc_it);
                    break;
                    
                case MAST::SURFACE_RADIATION_HEAT_FLUX:
                    surface_radiation_residual(request_jacobian,
                                               f, jac,
                                               it->first,
                                               **bc_it);
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





void
MAST::HeatConductionElementBase::
volume_external_residual (bool request_jacobian,
                          RealVectorX& f,
                          RealMatrixX& jac,
                          std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    typedef std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    
    libMesh::subdomain_id_type sid = _elem.get_reference_elem().subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::HEAT_FLUX:
                surface_flux_residual(request_jacobian,
                                      f, jac,
                                      *it.first->second);
                break;
            case MAST::CONVECTION_HEAT_FLUX:
                surface_convection_residual(request_jacobian,
                                            f, jac,
                                            *it.first->second);
                break;
                
            case MAST::SURFACE_RADIATION_HEAT_FLUX:
                surface_radiation_residual(request_jacobian,
                                           f, jac,
                                           *it.first->second);
                break;
                
            case MAST::HEAT_SOURCE:
                volume_heat_source_residual(request_jacobian,
                                            f, jac,
                                            *it.first->second);
                break;
                
            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    
}



void
MAST::HeatConductionElementBase::
side_external_residual_sensitivity (const MAST::FunctionBase& p,
                                    RealVectorX& f,
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
                case MAST::HEAT_FLUX:
                    surface_flux_residual_sensitivity(p,
                                                      f,
                                                      it->first,
                                                      **bc_it);
                    break;
                    
                case MAST::CONVECTION_HEAT_FLUX:
                    surface_convection_residual_sensitivity(p,
                                                            f,
                                                            it->first,
                                                            **bc_it);
                    break;
                    
                case MAST::SURFACE_RADIATION_HEAT_FLUX:
                    surface_radiation_residual_sensitivity(p,
                                                           f,
                                                           it->first,
                                                           **bc_it);
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




void
MAST::HeatConductionElementBase::
volume_external_residual_sensitivity (const MAST::FunctionBase& p,
                                      RealVectorX& f,
                                      std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    typedef std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    
    libMesh::subdomain_id_type sid = _elem.get_reference_elem().subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::HEAT_FLUX:
                surface_flux_residual_sensitivity(p,
                                                  f,
                                                  *it.first->second);
                break;

            case MAST::CONVECTION_HEAT_FLUX:
                surface_convection_residual_sensitivity(p,
                                                        f,
                                                        *it.first->second);
                break;
                
            case MAST::SURFACE_RADIATION_HEAT_FLUX:
                surface_radiation_residual_sensitivity(p,
                                                       f,
                                                       *it.first->second);
                break;
                
            case MAST::HEAT_SOURCE:
                volume_heat_source_residual_sensitivity(p,
                                                        f,
                                                        *it.first->second);
                break;
                
            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    
}



void
MAST::HeatConductionElementBase::
volume_external_residual_boundary_velocity(const MAST::FunctionBase& p,
                                           RealVectorX& f,
                                           const unsigned int s,
                                           const MAST::FieldFunction<RealVectorX>& vel_f,
                                           std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    typedef std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    
    libMesh::subdomain_id_type sid = _elem.get_reference_elem().subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::HEAT_FLUX:
                surface_flux_boundary_velocity(p,
                                               f,
                                               s,
                                               vel_f,
                                               *it.first->second);
                break;
                
            case MAST::CONVECTION_HEAT_FLUX:
                surface_convection_boundary_velocity(p,
                                                     f,
                                                     s,
                                                     vel_f,
                                                     *it.first->second);
                break;
                
            case MAST::SURFACE_RADIATION_HEAT_FLUX:
                surface_radiation_boundary_velocity(p,
                                                    f,
                                                    s,
                                                    vel_f,
                                                    *it.first->second);
                break;
                
            case MAST::HEAT_SOURCE:
                volume_heat_source_boundary_velocity(p,
                                                     f,
                                                     s,
                                                     vel_f,
                                                     *it.first->second);
                break;
                
            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    
}




void
MAST::HeatConductionElementBase::
internal_residual_sensitivity (const MAST::FunctionBase& p,
                               RealVectorX& f) {
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(true, false));

    const std::vector<Real>& JxW           = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = fe->get_xyz();
    const unsigned int
    n_phi  = fe->n_shape_functions(),
    dim    = _elem.dim();
    
    RealMatrixX
    material_mat   = RealMatrixX::Zero(dim, dim),
    dmaterial_mat  = RealMatrixX::Zero(dim, dim), // for calculation of Jac when k is temp. dep.
    mat_n2n2       = RealMatrixX::Zero(n_phi, n_phi);
    RealVectorX
    vec1     = RealVectorX::Zero(1),
    vec2_n2  = RealVectorX::Zero(n_phi),
    flux     = RealVectorX::Zero(dim);
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX> > conductance =
    _property.thermal_conductance_matrix(*this);
    
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix Bmat; // for calculation of Jac when k is temp. dep.
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_mass_fem_operator(qp, *fe, Bmat);
        Bmat.right_multiply(vec1, _sol);
        
        if (_active_sol_function)
            dynamic_cast<MAST::MeshFieldFunction*>
            (_active_sol_function)->set_element_quadrature_point_solution(vec1);
        
        conductance->derivative(p, xyz[qp], _time, material_mat);
        
        _initialize_fem_gradient_operator(qp, dim, *fe, dBmat);
        
        // calculate the flux for each dimension and add its weighted
        // component to the residual
        flux.setZero();
        for (unsigned int j=0; j<dim; j++) {
            dBmat[j].right_multiply(vec1, _sol);        // dT_dxj
            
            for (unsigned int i=0; i<dim; i++)
                flux(i) += vec1(0) * material_mat(i,j); // q_i = k_ij dT_dxj
        }
        
        // now add to the residual vector
        for (unsigned int i=0; i<dim; i++) {
            vec1(0)  = flux(i);
            dBmat[i].vector_mult_transpose(vec2_n2, vec1);
            f += JxW[qp] * vec2_n2;
        }
    }
    
    if (_active_sol_function)
        dynamic_cast<MAST::MeshFieldFunction*>
        (_active_sol_function)->clear_element_quadrature_point_solution();
    
}


void
MAST::HeatConductionElementBase::
internal_residual_boundary_velocity (const MAST::FunctionBase& p,
                                     RealVectorX& f,
                                     const unsigned int s,
                                     const MAST::FieldFunction<RealVectorX>& vel_f) {
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, true, false));

    std::vector<Real> JxW_Vn                        = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz          = fe->get_xyz();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals_for_local_coordinate();

    const unsigned int
    n_phi  = fe->n_shape_functions(),
    dim    = _elem.dim();
    
    RealMatrixX
    material_mat   = RealMatrixX::Zero(dim, dim),
    dmaterial_mat  = RealMatrixX::Zero(dim, dim), // for calculation of Jac when k is temp. dep.
    mat_n2n2       = RealMatrixX::Zero(n_phi, n_phi);
    RealVectorX
    vec1     = RealVectorX::Zero(1),
    vec2_n2  = RealVectorX::Zero(n_phi),
    flux     = RealVectorX::Zero(dim),
    vel      = RealVectorX::Zero(dim);
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX> > conductance =
    _property.thermal_conductance_matrix(*this);
    
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix Bmat; // for calculation of Jac when k is temp. dep.
    
    Real
    vn  = 0.;
    
    // modify the JxW_Vn by multiplying the normal velocity to it
    for (unsigned int qp=0; qp<JxW_Vn.size(); qp++) {
        
        vel_f(xyz[qp], _time, vel);
        vn = 0.;
        for (unsigned int i=0; i<dim; i++)
            vn += vel(i)*face_normals[qp](i);
        JxW_Vn[qp] *= vn;
    }

    
    for (unsigned int qp=0; qp<JxW_Vn.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_mass_fem_operator(qp, *fe, Bmat);
        Bmat.right_multiply(vec1, _sol);
        
        if (_active_sol_function)
            dynamic_cast<MAST::MeshFieldFunction*>
            (_active_sol_function)->set_element_quadrature_point_solution(vec1);
        
        (*conductance)(xyz[qp], _time, material_mat);
        
        _initialize_fem_gradient_operator(qp, dim, *fe, dBmat);
        
        // calculate the flux for each dimension and add its weighted
        // component to the residual
        flux.setZero();
        for (unsigned int j=0; j<dim; j++) {
            dBmat[j].right_multiply(vec1, _sol);        // dT_dxj
            
            for (unsigned int i=0; i<dim; i++)
                flux(i) += vec1(0) * material_mat(i,j); // q_i = k_ij dT_dxj
        }
        
        // now add to the residual vector
        for (unsigned int i=0; i<dim; i++) {
            vec1(0)  = flux(i);
            dBmat[i].vector_mult_transpose(vec2_n2, vec1);
            f += JxW_Vn[qp] * vec2_n2;
        }
    }
    
    if (_active_sol_function)
        dynamic_cast<MAST::MeshFieldFunction*>
        (_active_sol_function)->clear_element_quadrature_point_solution();
}



void
MAST::HeatConductionElementBase::
velocity_residual_sensitivity (const MAST::FunctionBase& p,
                               RealVectorX& f) {
    
    MAST::FEMOperatorMatrix Bmat;
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real>& JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz       = fe->get_xyz();
    
    const unsigned int
    n_phi      = fe->n_shape_functions(),
    dim        = _elem.dim();
    
    RealMatrixX
    material_mat    = RealMatrixX::Zero(dim, dim),
    mat_n2n2        = RealMatrixX::Zero(n_phi, n_phi),
    local_jac_xdot  = RealMatrixX::Zero(n_phi, n_phi);
    
    RealVectorX
    vec1    = RealVectorX::Zero(1),
    vec2_n2 = RealVectorX::Zero(n_phi),
    local_f = RealVectorX::Zero(n_phi);
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX> > capacitance =
    _property.thermal_capacitance_matrix(*this);
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        _initialize_mass_fem_operator(qp, *fe, Bmat);
        Bmat.right_multiply(vec1, _sol);               //  B * T
        
        if (_active_sol_function)
            dynamic_cast<MAST::MeshFieldFunction*>
            (_active_sol_function)->set_element_quadrature_point_solution(vec1);
        
        capacitance->derivative(p, xyz[qp], _time, material_mat);
        
        Bmat.right_multiply(vec1, _vel);               //  B * T_dot
        Bmat.vector_mult_transpose(vec2_n2, vec1);     //  B^T * B * T_dot
        
        local_f   += JxW[qp] * material_mat(0,0) * vec2_n2; // (rho*cp)*JxW B^T B T_dot
        
        if (_property.if_diagonal_mass_matrix()) {
            
            Bmat.right_multiply_transpose(mat_n2n2, Bmat);  // B^T B
            local_jac_xdot += JxW[qp] * material_mat(0,0) * mat_n2n2;  // B^T B * JxW (rho*cp)
        }
    }
    
    if (_property.if_diagonal_mass_matrix()) {
        
        for (unsigned int i=0; i<local_jac_xdot.rows(); i++) {
            
            Real a = local_jac_xdot.row(i).sum();
            local_jac_xdot.row(i).setZero();
            local_jac_xdot(i,i) = a;
        }

        f += local_jac_xdot * _vel;
    }
    else
        f += local_f;
    
    if (_active_sol_function)
        dynamic_cast<MAST::MeshFieldFunction*>
        (_active_sol_function)->clear_element_quadrature_point_solution();
    
}



void
MAST::HeatConductionElementBase::
velocity_residual_boundary_velocity (const MAST::FunctionBase& p,
                                     RealVectorX& f,
                                     const unsigned int s,
                                     const MAST::FieldFunction<RealVectorX>& vel_f) {
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    std::vector<Real> JxW_Vn                        = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz          = fe->get_xyz();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals_for_local_coordinate();

    MAST::FEMOperatorMatrix Bmat;
    
    const unsigned int
    n_phi      = fe->n_shape_functions(),
    dim        = _elem.dim();
    
    RealMatrixX
    material_mat    = RealMatrixX::Zero(dim, dim),
    mat_n2n2        = RealMatrixX::Zero(n_phi, n_phi),
    local_jac_xdot  = RealMatrixX::Zero(n_phi, n_phi);
    RealVectorX
    vec1    = RealVectorX::Zero(1),
    vec2_n2 = RealVectorX::Zero(n_phi),
    local_f = RealVectorX::Zero(n_phi),
    vel     = RealVectorX::Zero(dim);
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX> > capacitance =
    _property.thermal_capacitance_matrix(*this);

    Real
    vn  = 0.;
    
    // modify the JxW_Vn by multiplying the normal velocity to it
    for (unsigned int qp=0; qp<JxW_Vn.size(); qp++) {
        
        vel_f(xyz[qp], _time, vel);
        vn = 0.;
        for (unsigned int i=0; i<dim; i++)
            vn += vel(i)*face_normals[qp](i);
        JxW_Vn[qp] *= vn;
    }

    for (unsigned int qp=0; qp<JxW_Vn.size(); qp++) {
        
        _initialize_mass_fem_operator(qp, *fe, Bmat);
        Bmat.right_multiply(vec1, _sol);               //  B * T
        
        if (_active_sol_function)
            dynamic_cast<MAST::MeshFieldFunction*>
            (_active_sol_function)->set_element_quadrature_point_solution(vec1);
        
        (*capacitance)(xyz[qp], _time, material_mat);
        
        Bmat.right_multiply(vec1, _vel);               //  B * T_dot
        Bmat.vector_mult_transpose(vec2_n2, vec1);     //  B^T * B * T_dot
        
        local_f      += JxW_Vn[qp] * material_mat(0,0) * vec2_n2; // (rho*cp)*JxW B^T B T_dot

        if (_property.if_diagonal_mass_matrix()) {
            
            Bmat.right_multiply_transpose(mat_n2n2, Bmat);  // B^T B
            local_jac_xdot += JxW_Vn[qp] * material_mat(0,0) * mat_n2n2;  // B^T B * JxW (rho*cp)
        }
    }
    
    // diagonalize the matrix and compute the residual based on that.
    if (_property.if_diagonal_mass_matrix()) {

        for (unsigned int i=0; i<local_jac_xdot.rows(); i++) {
            
            Real a = local_jac_xdot.row(i).sum();
            local_jac_xdot.row(i).setZero();
            local_jac_xdot(i,i) = a;
        }

        f +=  local_jac_xdot * _vel;
    }
    else
        f += local_f;

    
    if (_active_sol_function)
        dynamic_cast<MAST::MeshFieldFunction*>
        (_active_sol_function)->clear_element_quadrature_point_solution();
}





void
MAST::HeatConductionElementBase::
surface_flux_residual(bool request_jacobian,
                      RealVectorX& f,
                      RealMatrixX& jac,
                      const unsigned int s,
                      MAST::BoundaryConditionBase& bc) {
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &func    = bc.get<MAST::FieldFunction<Real> >("heat_flux"),
    *section = _property.section(*this);
    
    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = fe->get_phi();
    const unsigned int n_phi                     = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    Real  flux, th;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get the value of flux = q_i . n_i
        func(qpoint[qp], _time, flux);
        if (section) (*section)(qpoint[qp], _time, th);
        else th  = 1.;

        f   +=  JxW[qp] * phi_vec * flux * th;
    }
}





void
MAST::HeatConductionElementBase::
surface_flux_residual(bool request_jacobian,
                      RealVectorX& f,
                      RealMatrixX& jac,
                      MAST::BoundaryConditionBase& bc) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>& func =
    bc.get<MAST::FieldFunction<Real> >("heat_flux");
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = fe->get_phi();
    const unsigned int n_phi                     = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    Real  flux;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get the value of flux = q_i . n_i
        func(qpoint[qp], _time, flux);
        
        f   +=  JxW[qp] * phi_vec * flux;
    }
}




void
MAST::HeatConductionElementBase::
surface_flux_residual_sensitivity(const MAST::FunctionBase& p,
                                  RealVectorX& f,
                                  const unsigned int s,
                                  MAST::BoundaryConditionBase& bc) {
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &func    = bc.get<MAST::FieldFunction<Real> >("heat_flux"),
    *section = _property.section(*this);

    
    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = fe->get_phi();
    const unsigned int n_phi                     = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    Real  flux, dflux, th, dth;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get the value of flux = q_i . n_i
        func(qpoint[qp], _time, flux);
        func.derivative(p, qpoint[qp], _time, dflux);
        if (section) {
            
            (*section)(qpoint[qp], _time, th);
            section->derivative(p, qpoint[qp], _time, dth);
        }
        else {
            th  = 1.;
            dth = 0.;
        }
        
        f   +=  JxW[qp] * phi_vec * (flux*dth + dflux*th);
    }
}




void
MAST::HeatConductionElementBase::
surface_flux_residual_sensitivity(const MAST::FunctionBase& p,
                                  RealVectorX& f,
                                  MAST::BoundaryConditionBase& bc) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>& func =
    bc.get<MAST::FieldFunction<Real> >("heat_flux");
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = fe->get_phi();
    const unsigned int n_phi                     = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    Real  flux;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get the value of flux = q_i . n_i
        func.derivative(p, qpoint[qp], _time, flux);
        
        f   +=  JxW[qp] * phi_vec * flux;
    }
}



void
MAST::HeatConductionElementBase::
surface_flux_boundary_velocity(const MAST::FunctionBase& p,
                               RealVectorX& f,
                               const unsigned int s,
                               const MAST::FieldFunction<RealVectorX>& vel_f,
                               MAST::BoundaryConditionBase& bc) {
    
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    std::vector<Real> JxW_Vn                        = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz          = fe->get_xyz();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals_for_local_coordinate();
    const std::vector<std::vector<Real> >& phi      = fe->get_phi();
    const unsigned int
    n_phi    = (unsigned int)phi.size(),
    dim      = _elem.dim();

    RealVectorX
    phi_vec  = RealVectorX::Zero(n_phi),
    vel      = RealVectorX::Zero(dim);
    Real  flux;

    Real
    vn  = 0.;
    
    // modify the JxW_Vn by multiplying the normal velocity to it
    for (unsigned int qp=0; qp<JxW_Vn.size(); qp++) {
        
        vel_f(xyz[qp], _time, vel);
        vn = 0.;
        for (unsigned int i=0; i<dim; i++)
            vn += vel(i)*face_normals[qp](i);
        JxW_Vn[qp] *= vn;
    }

    // get the function from this boundary condition
    const MAST::FieldFunction<Real>& func =
    bc.get<MAST::FieldFunction<Real> >("heat_flux");
    

    for (unsigned int qp=0; qp<xyz.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get the value of flux = q_i . n_i
        func(xyz[qp], _time, flux);
        
        f   +=  JxW_Vn[qp] * phi_vec * flux;
    }
}


void
MAST::HeatConductionElementBase::
surface_convection_residual(bool request_jacobian,
                            RealVectorX& f,
                            RealMatrixX& jac,
                            const unsigned int s,
                            MAST::BoundaryConditionBase& bc) {
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &coeff   = bc.get<MAST::FieldFunction<Real> >("convection_coeff"),
    &T_amb   = bc.get<MAST::FieldFunction<Real> >("ambient_temperature"),
    *section = _property.section(*this);

    const std::vector<Real> &JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int n_phi                   = (unsigned int)phi.size();
    
    
    RealVectorX  phi_vec  = RealVectorX::Zero(n_phi);
    RealMatrixX  mat      = RealMatrixX::Zero(n_phi, n_phi);
    Real temp, amb_temp, h_coeff, th;
    MAST::FEMOperatorMatrix Bmat;
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        coeff(qpoint[qp], _time, h_coeff);
        T_amb(qpoint[qp], _time, amb_temp);
        temp  = phi_vec.dot(_sol);
        if (section) (*section)(qpoint[qp], _time, th);
        else th  = 1.;

        // normal flux is given as:
        // qi_ni = h_coeff * (T - T_amb)
        //
        f   += JxW[qp] * phi_vec * th* h_coeff * (temp - amb_temp);
        
        if (request_jacobian) {
            
            Bmat.reinit(1, phi_vec);
            Bmat.right_multiply_transpose(mat, Bmat);
            jac += JxW[qp] * mat * th * h_coeff;
        }
    }
    
}





void
MAST::HeatConductionElementBase::
surface_convection_residual(bool request_jacobian,
                            RealVectorX& f,
                            RealMatrixX& jac,
                            MAST::BoundaryConditionBase& bc) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &coeff = bc.get<MAST::FieldFunction<Real> >("convection_coeff"),
    &T_amb = bc.get<MAST::FieldFunction<Real> >("ambient_temperature");
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real> &JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int n_phi                   = (unsigned int)phi.size();
    
    
    RealVectorX  phi_vec  = RealVectorX::Zero(n_phi);
    RealMatrixX  mat      = RealMatrixX::Zero(n_phi, n_phi);
    Real temp, amb_temp, h_coeff;
    MAST::FEMOperatorMatrix Bmat;
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        coeff(qpoint[qp], _time, h_coeff);
        T_amb(qpoint[qp], _time, amb_temp);
        temp  = phi_vec.dot(_sol);
        
        // normal flux is given as:
        // qi_ni = h_coeff * (T - T_amb)
        //
        f   += JxW[qp] * phi_vec * h_coeff * (temp - amb_temp);
        
        if (request_jacobian) {
            
            Bmat.reinit(1, phi_vec);
            Bmat.right_multiply_transpose(mat, Bmat);
            jac += JxW[qp] * mat * h_coeff;
        }
    }
    
}





void
MAST::HeatConductionElementBase::
surface_convection_residual_sensitivity(const MAST::FunctionBase& p,
                                        RealVectorX& f,
                                        const unsigned int s,
                                        MAST::BoundaryConditionBase& bc) {
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &coeff   = bc.get<MAST::FieldFunction<Real> >("convection_coeff"),
    &T_amb   = bc.get<MAST::FieldFunction<Real> >("ambient_temperature"),
    *section = _property.section(*this);
    
    const std::vector<Real> &JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int n_phi                   = (unsigned int)phi.size();
    
    
    RealVectorX  phi_vec  = RealVectorX::Zero(n_phi);
    RealMatrixX  mat      = RealMatrixX::Zero(n_phi, n_phi);
    Real
    temp, amb_temp, h_coeff, th,
    damb_temp, dh_coeff, dth;
    MAST::FEMOperatorMatrix Bmat;
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        coeff(qpoint[qp], _time, h_coeff);
        T_amb(qpoint[qp], _time, amb_temp);
        coeff.derivative(p, qpoint[qp], _time, dh_coeff);
        T_amb.derivative(p, qpoint[qp], _time, damb_temp);
        if (section) {
            
            (*section)(qpoint[qp], _time, th);
            section->derivative(p, qpoint[qp], _time, dth);
        }
        else {
            th  = 1.;
            dth = 0.;
        }
        temp  = phi_vec.dot(_sol);
        
        // normal flux is given as:
        // qi_ni = h_coeff * (T - T_amb)
        //
        f   += JxW[qp] * phi_vec * (th * dh_coeff * (temp - amb_temp) +
                                    th * h_coeff * (-damb_temp) +
                                    dth* h_coeff * (temp - amb_temp));
    }
}




void
MAST::HeatConductionElementBase::
surface_convection_residual_sensitivity(const MAST::FunctionBase& p,
                                        RealVectorX& f,
                                        MAST::BoundaryConditionBase& bc) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &coeff = bc.get<MAST::FieldFunction<Real> >("convection_coeff"),
    &T_amb = bc.get<MAST::FieldFunction<Real> >("ambient_temperature");

    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real> &JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int n_phi                   = (unsigned int)phi.size();
    
    
    RealVectorX  phi_vec  = RealVectorX::Zero(n_phi);
    RealMatrixX  mat      = RealMatrixX::Zero(n_phi, n_phi);
    Real
    temp, amb_temp, h_coeff,
    damb_temp, dh_coeff;
    MAST::FEMOperatorMatrix Bmat;
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        coeff(qpoint[qp], _time, h_coeff);
        T_amb(qpoint[qp], _time, amb_temp);
        coeff.derivative(p, qpoint[qp], _time, dh_coeff);
        T_amb.derivative(p, qpoint[qp], _time, damb_temp);
        temp  = phi_vec.dot(_sol);
        
        // normal flux is given as:
        // qi_ni = h_coeff * (T - T_amb)
        //
        f   += JxW[qp] * phi_vec * (dh_coeff * (temp - amb_temp)+
                                    h_coeff * (-damb_temp));
    }
    
}



void
MAST::HeatConductionElementBase::
surface_convection_boundary_velocity(const MAST::FunctionBase& p,
                                     RealVectorX& f,
                                     const unsigned int s,
                                     const MAST::FieldFunction<RealVectorX>& vel_f,
                                     MAST::BoundaryConditionBase& bc) {
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    std::vector<Real> JxW_Vn                        = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz          = fe->get_xyz();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals_for_local_coordinate();
    const std::vector<std::vector<Real> >& phi      = fe->get_phi();

    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &coeff = bc.get<MAST::FieldFunction<Real> >("convection_coeff"),
    &T_amb = bc.get<MAST::FieldFunction<Real> >("ambient_temperature");
    
    const unsigned int
    n_phi   = (unsigned int)phi.size(),
    dim     = _elem.dim();
    
    
    RealVectorX
    phi_vec  = RealVectorX::Zero(n_phi),
    vel      = RealVectorX::Zero(dim);
    RealMatrixX
    mat      = RealMatrixX::Zero(n_phi, n_phi);
    
    Real temp, amb_temp, h_coeff;
    MAST::FEMOperatorMatrix Bmat;
    
    Real
    vn  = 0.;
    
    
    // modify the JxW_Vn by multiplying the normal velocity to it
    for (unsigned int qp=0; qp<JxW_Vn.size(); qp++) {
        
        vel_f(xyz[qp], _time, vel);
        vn = 0.;
        for (unsigned int i=0; i<dim; i++)
            vn += vel(i)*face_normals[qp](i);
        JxW_Vn[qp] *= vn;
    }

    for (unsigned int qp=0; qp<xyz.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        coeff(xyz[qp], _time, h_coeff);
        T_amb(xyz[qp], _time, amb_temp);
        temp  = phi_vec.dot(_sol);
        
        // normal flux is given as:
        // qi_ni = h_coeff * (T - T_amb)
        //
        f   += JxW_Vn[qp] * phi_vec * h_coeff * (temp - amb_temp);
    }
    
}




void
MAST::HeatConductionElementBase::
surface_radiation_residual(bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac,
                           const unsigned int s,
                           MAST::BoundaryConditionBase& bc) {
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &emissivity = bc.get<MAST::FieldFunction<Real> >("emissivity"),
    *section    = _property.section(*this);
    
    const MAST::Parameter
    &T_amb      = bc.get<MAST::Parameter>("ambient_temperature"),
    &T_ref_zero = bc.get<MAST::Parameter>("reference_zero_temperature"),
    &sb_const   = bc.get<MAST::Parameter>("stefan_bolzmann_constant");
    
    
    const std::vector<Real> &JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int n_phi                   = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    RealMatrixX mat      = RealMatrixX::Zero(n_phi, n_phi);
    const Real
    sbc      = sb_const(),
    amb_temp = T_amb(),
    zero_ref = T_ref_zero();
    Real temp, emiss, th;
    MAST::FEMOperatorMatrix Bmat;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        emissivity(qpoint[qp], _time, emiss);
        temp  = phi_vec.dot(_sol);
        if (section) (*section)(qpoint[qp], _time, th);
        else th  = 1.;

        f   += JxW[qp] * phi_vec * sbc * emiss * th *
        (pow(temp+zero_ref, 4.) - pow(amb_temp+zero_ref, 4.));
        
        if (request_jacobian) {
            
            Bmat.reinit(1, phi_vec);
            Bmat.right_multiply_transpose(mat, Bmat);
            jac +=  JxW[qp] * mat * sbc * emiss * th * 4. * pow(temp+zero_ref, 3.);
        }
    }
}




void
MAST::HeatConductionElementBase::
surface_radiation_residual(bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac,
                           MAST::BoundaryConditionBase& bc) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &emissivity = bc.get<MAST::FieldFunction<Real> >("emissivity");
    
    const MAST::Parameter
    &T_amb      = bc.get<MAST::Parameter>("ambient_temperature"),
    &T_ref_zero = bc.get<MAST::Parameter>("reference_zero_temperature"),
    &sb_const   = bc.get<MAST::Parameter>("stefan_bolzmann_constant");
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real> &JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int n_phi                   = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    RealMatrixX mat      = RealMatrixX::Zero(n_phi, n_phi);
    const Real
    sbc      = sb_const(),
    amb_temp = T_amb(),
    zero_ref = T_ref_zero();
    Real temp, emiss;
    MAST::FEMOperatorMatrix Bmat;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        emissivity(qpoint[qp], _time, emiss);
        temp  = phi_vec.dot(_sol);
        
        f   += JxW[qp] * phi_vec * sbc * emiss *
        (pow(temp+zero_ref, 4.) - pow(amb_temp+zero_ref, 4.));
        
        if (request_jacobian) {
            
            Bmat.reinit(1, phi_vec);
            Bmat.right_multiply_transpose(mat, Bmat);
            jac +=  JxW[qp] * mat * sbc * emiss * 4. * pow(temp+zero_ref, 3.);
        }
    }
}





void
MAST::HeatConductionElementBase::
surface_radiation_residual_sensitivity(const MAST::FunctionBase& p,
                                       RealVectorX& f,
                                       const unsigned int s,
                                       MAST::BoundaryConditionBase& bc) {

    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &emissivity = bc.get<MAST::FieldFunction<Real> >("emissivity"),
    *section    = _property.section(*this);
    
    const MAST::Parameter
    &T_amb      = bc.get<MAST::Parameter>("ambient_temperature"),
    &T_ref_zero = bc.get<MAST::Parameter>("reference_zero_temperature"),
    &sb_const   = bc.get<MAST::Parameter>("stefan_bolzmann_constant");
    
    
    const std::vector<Real> &JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int n_phi                   = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    RealMatrixX mat      = RealMatrixX::Zero(n_phi, n_phi);
    const Real
    sbc      = sb_const(),
    amb_temp = T_amb(),
    zero_ref = T_ref_zero();
    Real temp, emiss, demiss, th, dth;
    MAST::FEMOperatorMatrix Bmat;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        emissivity(qpoint[qp], _time, emiss);
        emissivity.derivative(p, qpoint[qp], _time, demiss);
        temp  = phi_vec.dot(_sol);
        if (section) {
            
            (*section)(qpoint[qp], _time, th);
            section->derivative(p, qpoint[qp], _time, dth);
        }
        else {
            th  = 1.;
            dth = 0.;
        }

        f   += JxW[qp] * phi_vec * sbc * (demiss * th + emiss * dth) *
        (pow(temp+zero_ref, 4.) - pow(amb_temp+zero_ref, 4.));
    }
}




void
MAST::HeatConductionElementBase::
surface_radiation_residual_sensitivity(const MAST::FunctionBase& p,
                                       RealVectorX& f,
                                       MAST::BoundaryConditionBase& bc) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &emissivity = bc.get<MAST::FieldFunction<Real> >("emissivity");
    
    const MAST::Parameter
    &T_amb      = bc.get<MAST::Parameter>("ambient_temperature"),
    &T_ref_zero = bc.get<MAST::Parameter>("reference_zero_temperature"),
    &sb_const   = bc.get<MAST::Parameter>("stefan_bolzmann_constant");
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real> &JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int n_phi                   = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    RealMatrixX mat      = RealMatrixX::Zero(n_phi, n_phi);
    const Real
    sbc      = sb_const(),
    amb_temp = T_amb(),
    zero_ref = T_ref_zero();
    Real temp, demiss;
    MAST::FEMOperatorMatrix Bmat;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        emissivity.derivative(p, qpoint[qp], _time, demiss);
        temp  = phi_vec.dot(_sol);
        
        f   += JxW[qp] * phi_vec * sbc * demiss *
        (pow(temp+zero_ref, 4.) - pow(amb_temp+zero_ref, 4.));
    }
}


void
MAST::HeatConductionElementBase::
surface_radiation_boundary_velocity(const MAST::FunctionBase& p,
                                    RealVectorX& f,
                                    const unsigned int s,
                                    const MAST::FieldFunction<RealVectorX>& vel_f,
                                    MAST::BoundaryConditionBase& bc) {
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    std::vector<Real> JxW_Vn                        = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz          = fe->get_xyz();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals_for_local_coordinate();
    const std::vector<std::vector<Real> >& phi      = fe->get_phi();

    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &emissivity = bc.get<MAST::FieldFunction<Real> >("emissivity");
    
    const MAST::Parameter
    &T_amb      = bc.get<MAST::Parameter>("ambient_temperature"),
    &T_ref_zero = bc.get<MAST::Parameter>("reference_zero_temperature"),
    &sb_const   = bc.get<MAST::Parameter>("stefan_bolzmann_constant");
    
    const unsigned int
    n_phi       = (unsigned int)phi.size(),
    dim         = _elem.dim();
    
    RealVectorX
    phi_vec  = RealVectorX::Zero(n_phi),
    vel      = RealVectorX::Zero(dim);
    RealMatrixX
    mat      = RealMatrixX::Zero(n_phi, n_phi);
    
    const Real
    sbc      = sb_const(),
    amb_temp = T_amb(),
    zero_ref = T_ref_zero();
    
    Real
    temp   = 0.,
    emiss  = 0.,
    vn     = 0.;
    
    MAST::FEMOperatorMatrix Bmat;
    
    
    // modify the JxW_Vn by multiplying the normal velocity to it
    for (unsigned int qp=0; qp<JxW_Vn.size(); qp++) {
        
        vel_f(xyz[qp], _time, vel);
        vn = 0.;
        for (unsigned int i=0; i<dim; i++)
            vn += vel(i)*face_normals[qp](i);
        JxW_Vn[qp] *= vn;
    }

    for (unsigned int qp=0; qp<xyz.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        emissivity(xyz[qp], _time, emiss);
        temp  = phi_vec.dot(_sol);
        
        f   += JxW_Vn[qp] * phi_vec * sbc * emiss *
        (pow(temp+zero_ref, 4.) - pow(amb_temp+zero_ref, 4.));
    }
}



void
MAST::HeatConductionElementBase::
volume_heat_source_residual(bool request_jacobian,
                            RealVectorX& f,
                            RealMatrixX& jac,
                            MAST::BoundaryConditionBase& bc) {
    
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &func       = bc.get<MAST::FieldFunction<Real> >("heat_source"),
    *section    = _property.section(*this);
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));
    
    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = fe->get_phi();
    const unsigned int n_phi                     = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    Real  source, th;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get the value of heat source
        func(qpoint[qp], _time, source);
        if (section) (*section)(qpoint[qp], _time, th);
        else th  = 1.;

        f   -=  JxW[qp] * phi_vec * source * th;
    }
}



void
MAST::HeatConductionElementBase::
volume_heat_source_residual_sensitivity(const MAST::FunctionBase& p,
                                        RealVectorX& f,
                                        MAST::BoundaryConditionBase& bc) {

    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &func       = bc.get<MAST::FieldFunction<Real> >("heat_source"),
    *section    = _property.section(*this);
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(false, false));

    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = fe->get_phi();
    const unsigned int n_phi                     = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    Real  source, dsource, th, dth;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get the value of heat source
        func(qpoint[qp], _time, source);
        func.derivative(p, qpoint[qp], _time, dsource);
        if (section) {
            
            (*section)(qpoint[qp], _time, th);
            section->derivative(p, qpoint[qp], _time, dth);
        }
        else {
            th  = 1.;
            dth = 0.;
        }

        f   -=  JxW[qp] * phi_vec * (dsource * th + source * dth);
    }
}


void
MAST::HeatConductionElementBase::
volume_heat_source_boundary_velocity(const MAST::FunctionBase& p,
                                     RealVectorX& f,
                                     const unsigned int s,
                                     const MAST::FieldFunction<RealVectorX>& vel_f,
                                     MAST::BoundaryConditionBase& bc) {
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, false, false));

    std::vector<Real> JxW_Vn                        = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz          = fe->get_xyz();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals_for_local_coordinate();
    const std::vector<std::vector<Real> >& phi      = fe->get_phi();

    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &func       = bc.get<MAST::FieldFunction<Real> >("heat_source"),
    *section    = _property.section(*this);
    
    const unsigned int
    n_phi   = (unsigned int)phi.size(),
    dim     = _elem.dim();
    
    RealVectorX
    phi_vec  = RealVectorX::Zero(n_phi),
    vel      = RealVectorX::Zero(dim);
    Real
    source = 0.,
    th     = 0.,
    vn     = 0.;
    
    // modify the JxW_Vn by multiplying the normal velocity to it
    for (unsigned int qp=0; qp<JxW_Vn.size(); qp++) {
        
        vel_f(xyz[qp], _time, vel);
        vn = 0.;
        for (unsigned int i=0; i<dim; i++)
            vn += vel(i)*face_normals[qp](i);
        JxW_Vn[qp] *= vn;
    }

    for (unsigned int qp=0; qp<xyz.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get the value of heat source
        func(xyz[qp], _time, source);
        if (section) (*section)(xyz[qp], _time, th);
        else th  = 1.;

        f   -=  JxW_Vn[qp] * phi_vec * source * th;
    }
}




void
MAST::HeatConductionElementBase::
_initialize_mass_fem_operator(const unsigned int qp,
                              const MAST::FEBase& fe,
                              MAST::FEMOperatorMatrix& Bmat) {
    
    const std::vector<std::vector<Real> >& phi_fe = fe.get_phi();
    
    const unsigned int n_phi = (unsigned int)phi_fe.size();
    
    RealVectorX phi = RealVectorX::Zero(n_phi);
    
    // shape function values
    // N
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = phi_fe[i_nd][qp];
    
    Bmat.reinit(1, phi);
}




void
MAST::HeatConductionElementBase::
_initialize_fem_gradient_operator(const unsigned int qp,
                                  const unsigned int dim,
                                  const MAST::FEBase& fe,
                                  std::vector<MAST::FEMOperatorMatrix>& dBmat) {
    
    libmesh_assert(dBmat.size() == dim);
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe.get_dphi();
    
    const unsigned int n_phi = (unsigned int)dphi.size();
    RealVectorX phi = RealVectorX::Zero(n_phi);
    
    // now set the shape function values
    for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
        
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi(i_nd) = dphi[i_nd][qp](i_dim);
        dBmat[i_dim].reinit(1, phi); //  dT/dx_i
    }
}

