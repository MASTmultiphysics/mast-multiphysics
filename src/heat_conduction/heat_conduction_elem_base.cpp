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
#include "heat_conduction/heat_conduction_elem_base.h"
#include "numerics/fem_operator_matrix.h"
#include "base/system_initialization.h"
#include "base/field_function_base.h"
#include "base/parameter.h"
#include "base/boundary_condition_base.h"
#include "property_cards/element_property_card_base.h"
#include "property_cards/element_property_card_1D.h"
#include "mesh/local_elem_base.h"
#include "mesh/local_1d_elem.h"
#include "mesh/local_2d_elem.h"
#include "mesh/local_3d_elem.h"
#include "base/mesh_field_function.h"


MAST::HeatConductionElementBase::
HeatConductionElementBase(MAST::SystemInitialization& sys,
                          const libMesh::Elem& elem,
                          const MAST::ElementPropertyCardBase& p):
MAST::ElementBase(sys, elem),
_property(p) {

    MAST::LocalElemBase* rval = nullptr;
    
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
    
    // now initialize the finite element data structures
    _init_fe_and_qrule(get_elem_for_quadrature(), &_fe, &_qrule);
}



MAST::HeatConductionElementBase::~HeatConductionElementBase() {
    
}




bool
MAST::HeatConductionElementBase::internal_residual (bool request_jacobian,
                                                    RealVectorX& f,
                                                    RealMatrixX& jac) {
    
    const std::vector<Real>& JxW           = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int
    n_phi  = _fe->n_shape_functions(),
    dim    = _elem.dim();
    
    RealMatrixX
    material_mat   = RealMatrixX::Zero(dim, dim),
    dmaterial_mat  = RealMatrixX::Zero(dim, dim), // for calculation of Jac when k is temp. dep.
    mat_n2n2       = RealMatrixX::Zero(n_phi, n_phi);
    RealVectorX
    vec1     = RealVectorX::Zero(1),
    vec2_n2  = RealVectorX::Zero(n_phi),
    flux     = RealVectorX::Zero(dim);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > conductance =
    _property.thermal_conductance_matrix(*this);
    
    libMesh::Point p;
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix Bmat; // for calculation of Jac when k is temp. dep.

    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {

        // initialize the Bmat operator for this term
        _initialize_mass_fem_operator(qp, *_fe, Bmat);
        Bmat.right_multiply(vec1, _sol);

        if (_active_sol_function)
            dynamic_cast<MAST::MeshFieldFunction<RealVectorX>*>
            (_active_sol_function)->set_element_quadrature_point_solution(vec1);

        _local_elem->global_coordinates_location(xyz[qp], p);
        
        (*conductance)(p, _time, material_mat);

        _initialize_fem_gradient_operator(qp, dim, *_fe, dBmat);
        
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
                conductance->derivative(MAST::PARTIAL_DERIVATIVE,
                                        *_active_sol_function,
                                        p,
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
        dynamic_cast<MAST::MeshFieldFunction<RealVectorX>*>
        (_active_sol_function)->clear_element_quadrature_point_solution();

    return request_jacobian;
}





bool
MAST::HeatConductionElementBase::velocity_residual (bool request_jacobian,
                                                    RealVectorX& f,
                                                    RealMatrixX& jac_xdot,
                                                    RealMatrixX& jac) {
    MAST::FEMOperatorMatrix Bmat;
    
    const std::vector<Real>& JxW                 = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz       = _fe->get_xyz();
    
    const unsigned int
    n_phi      = _fe->n_shape_functions(),
    dim        = _elem.dim();
    
    RealMatrixX
    material_mat    = RealMatrixX::Zero(dim, dim),
    mat_n2n2        = RealMatrixX::Zero(n_phi, n_phi);
    RealVectorX
    vec1    = RealVectorX::Zero(1),
    vec2_n2 = RealVectorX::Zero(n_phi);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > capacitance =
    _property.thermal_capacitance_matrix(*this);
    
    libMesh::Point p;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        _initialize_mass_fem_operator(qp, *_fe, Bmat);
        Bmat.right_multiply(vec1, _sol);               //  B * T
        
        if (_active_sol_function)
            dynamic_cast<MAST::MeshFieldFunction<RealVectorX>*>
            (_active_sol_function)->set_element_quadrature_point_solution(vec1);

        _local_elem->global_coordinates_location(xyz[qp], p);
        
        (*capacitance)(p, _time, material_mat);
        
        Bmat.right_multiply(vec1, _vel);               //  B * T_dot
        Bmat.vector_mult_transpose(vec2_n2, vec1);     //  B^T * B * T_dot
        
        f      += JxW[qp] * material_mat(0,0) * vec2_n2; // (rho*cp)*JxW B^T B T_dot
        
        if (request_jacobian) {
            
            Bmat.right_multiply_transpose(mat_n2n2, Bmat);  // B^T B
            jac_xdot += JxW[qp] * material_mat(0,0) * mat_n2n2;  // B^T B * JxW (rho*cp)
            
            // Jacobian contribution from int_omega B T d(rho*cp)/dT B
            if (_active_sol_function) {
                // get derivative of the conductance matrix wrt temperature
                capacitance->derivative(MAST::PARTIAL_DERIVATIVE,
                                        *_active_sol_function,
                                        p,
                                        _time, material_mat);
                
                if (material_mat(0,0) != 0.) { // no need to process for zero terms
                    
                    // B^T (T d(rho cp)/dT) B
                    jac += JxW[qp] * vec1(0) * material_mat(0,0) * mat_n2n2;
                }
            }
        }
    }
    
    
    if (_active_sol_function)
        dynamic_cast<MAST::MeshFieldFunction<RealVectorX>*>
        (_active_sol_function)->clear_element_quadrature_point_solution();

    return request_jacobian;
}




bool
MAST::HeatConductionElementBase::
side_external_residual (bool request_jacobian,
                        RealVectorX& f,
                        RealMatrixX& jac,
                        std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    typedef std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    const libMesh::BoundaryInfo& binfo = *_system.system().get_mesh().boundary_info;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    
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
                    case MAST::HEAT_FLUX:
                        surface_flux_residual(request_jacobian,
                                              f, jac,
                                              n,
                                              *it.first->second);
                        break;
                        
                    case MAST::CONVECTION_HEAT_FLUX:
                        surface_convection_residual(request_jacobian,
                                                    f, jac,
                                                    n,
                                                    *it.first->second);
                        break;
                        
                    case MAST::SURFACE_RADIATION_HEAT_FLUX:
                        surface_radiation_residual(request_jacobian,
                                                   f, jac,
                                                   n,
                                                   *it.first->second);
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
    return request_jacobian;
}





bool
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
    
    libMesh::subdomain_id_type sid = _elem.subdomain_id();
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
    
    return request_jacobian;
}



bool
MAST::HeatConductionElementBase::
side_external_residual_sensitivity (bool request_jacobian,
                                    RealVectorX& f,
                                    RealMatrixX& jac,
                                    std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    
    typedef std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    const libMesh::BoundaryInfo& binfo = *_system.system().get_mesh().boundary_info;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    
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
                    case MAST::HEAT_FLUX:
                        surface_flux_residual_sensitivity(request_jacobian,
                                                          f, jac,
                                                          n,
                                                          *it.first->second);
                        break;
                        
                    case MAST::CONVECTION_HEAT_FLUX:
                        surface_convection_residual_sensitivity(request_jacobian,
                                                                f, jac,
                                                                n,
                                                                *it.first->second);
                        break;
                        
                    case MAST::SURFACE_RADIATION_HEAT_FLUX:
                        surface_radiation_residual_sensitivity(request_jacobian,
                                                               f, jac,
                                                               n,
                                                               *it.first->second);
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
    return request_jacobian;
}




bool
MAST::HeatConductionElementBase::
volume_external_residual_sensitivity (bool request_jacobian,
                                      RealVectorX& f,
                                      RealMatrixX& jac,
                                      std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    typedef std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    
    libMesh::subdomain_id_type sid = _elem.subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {
                
            case MAST::HEAT_FLUX:
                surface_flux_residual_sensitivity(request_jacobian,
                                                  f, jac,
                                                  *it.first->second);
                break;

            case MAST::CONVECTION_HEAT_FLUX:
                surface_convection_residual_sensitivity(request_jacobian,
                                                        f, jac,
                                                        *it.first->second);
                break;
                
            case MAST::SURFACE_RADIATION_HEAT_FLUX:
                surface_radiation_residual_sensitivity(request_jacobian,
                                                       f, jac,
                                                       *it.first->second);
                break;
                
            case MAST::HEAT_SOURCE:
                volume_heat_source_residual_sensitivity(request_jacobian,
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
MAST::HeatConductionElementBase::
internal_residual_sensitivity (bool request_jacobian,
                               RealVectorX& f,
                               RealMatrixX& jac) {
    
    return request_jacobian;
}



bool
MAST::HeatConductionElementBase::
velocity_residual_sensitivity (bool request_jacobian,
                               RealVectorX& f,
                               RealMatrixX& jac) {
    
    return request_jacobian;
}








bool
MAST::HeatConductionElementBase::
surface_flux_residual(bool request_jacobian,
                      RealVectorX& f,
                      RealMatrixX& jac,
                      const unsigned int s,
                      MAST::BoundaryConditionBase& p) {
    
    // prepare the side finite element
    libMesh::FEBase *fe_ptr    = nullptr;
    libMesh::QBase  *qrule_ptr = nullptr;
    _get_side_fe_and_qrule(get_elem_for_quadrature(), s, &fe_ptr, &qrule_ptr, false);
    std::auto_ptr<libMesh::FEBase> fe(fe_ptr);
    std::auto_ptr<libMesh::QBase>  qrule(qrule_ptr);
    
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>& func =
    p.get<MAST::FieldFunction<Real> >("heat_flux");
    
    
    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = fe->get_phi();
    const unsigned int n_phi                     = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    libMesh::Point pt;
    Real  flux;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location (qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get the value of flux = q_i . n_i
        func(pt, _time, flux);
        
        f   +=  JxW[qp] * phi_vec * flux;
    }
    
    // calculation of the load vector is independent of solution
    return false;
}





bool
MAST::HeatConductionElementBase::
surface_flux_residual(bool request_jacobian,
                      RealVectorX& f,
                      RealMatrixX& jac,
                      MAST::BoundaryConditionBase& p) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>& func =
    p.get<MAST::FieldFunction<Real> >("heat_flux");
    
    
    const std::vector<Real> &JxW                 = _fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = _fe->get_phi();
    const unsigned int n_phi                     = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    libMesh::Point pt;
    Real  flux;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location (qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get the value of flux = q_i . n_i
        func(pt, _time, flux);
        
        f   +=  JxW[qp] * phi_vec * flux;
    }
    
    // calculation of the load vector is independent of solution
    return false;
}




bool
MAST::HeatConductionElementBase::
surface_flux_residual_sensitivity(bool request_jacobian,
                                  RealVectorX& f,
                                  RealMatrixX& jac,
                                  const unsigned int s,
                                  MAST::BoundaryConditionBase& p) {
    
    return false;
}




bool
MAST::HeatConductionElementBase::
surface_flux_residual_sensitivity(bool request_jacobian,
                                  RealVectorX& f,
                                  RealMatrixX& jac,
                                  MAST::BoundaryConditionBase& p) {
    
    return false;
}




bool
MAST::HeatConductionElementBase::
surface_convection_residual(bool request_jacobian,
                            RealVectorX& f,
                            RealMatrixX& jac,
                            const unsigned int s,
                            MAST::BoundaryConditionBase& p) {
    
    // prepare the side finite element
    libMesh::FEBase *fe_ptr    = nullptr;
    libMesh::QBase  *qrule_ptr = nullptr;
    _get_side_fe_and_qrule(get_elem_for_quadrature(), s, &fe_ptr, &qrule_ptr, false);
    std::auto_ptr<libMesh::FEBase> fe(fe_ptr);
    std::auto_ptr<libMesh::QBase>  qrule(qrule_ptr);

    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &coeff = p.get<MAST::FieldFunction<Real> >("convection_coeff"),
    &T_amb = p.get<MAST::FieldFunction<Real> >("ambient_temperature");
    
    const std::vector<Real> &JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const unsigned int n_phi                   = (unsigned int)phi.size();
    
    
    RealVectorX  phi_vec  = RealVectorX::Zero(n_phi);
    RealMatrixX  mat      = RealMatrixX::Zero(n_phi, n_phi);
    Real temp, amb_temp, h_coeff;
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location (qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        coeff(pt, _time, h_coeff);
        T_amb(pt, _time, amb_temp);
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
    
    return request_jacobian;
}





bool
MAST::HeatConductionElementBase::
surface_convection_residual(bool request_jacobian,
                            RealVectorX& f,
                            RealMatrixX& jac,
                            MAST::BoundaryConditionBase& p) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &coeff = p.get<MAST::FieldFunction<Real> >("convection_coeff"),
    &T_amb = p.get<MAST::FieldFunction<Real> >("ambient_temperature");
    
    const std::vector<Real> &JxW               = _fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = _fe->get_phi();
    const unsigned int n_phi                   = (unsigned int)phi.size();
    
    
    RealVectorX  phi_vec  = RealVectorX::Zero(n_phi);
    RealMatrixX  mat      = RealMatrixX::Zero(n_phi, n_phi);
    Real temp, amb_temp, h_coeff;
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location (qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        coeff(pt, _time, h_coeff);
        T_amb(pt, _time, amb_temp);
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
    
    return request_jacobian;
}





bool
MAST::HeatConductionElementBase::
surface_convection_residual_sensitivity(bool request_jacobian,
                                        RealVectorX& f,
                                        RealMatrixX& jac,
                                        const unsigned int s,
                                        MAST::BoundaryConditionBase& p) {
    
    return false;
}




bool
MAST::HeatConductionElementBase::
surface_convection_residual_sensitivity(bool request_jacobian,
                                        RealVectorX& f,
                                        RealMatrixX& jac,
                                        MAST::BoundaryConditionBase& p) {
    
    return false;
}





bool
MAST::HeatConductionElementBase::
surface_radiation_residual(bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac,
                           const unsigned int s,
                           MAST::BoundaryConditionBase& p) {
    
    // prepare the side finite element
    libMesh::FEBase *fe_ptr    = nullptr;
    libMesh::QBase  *qrule_ptr = nullptr;
    _get_side_fe_and_qrule(get_elem_for_quadrature(), s, &fe_ptr, &qrule_ptr, false);
    std::auto_ptr<libMesh::FEBase> fe(fe_ptr);
    std::auto_ptr<libMesh::QBase>  qrule(qrule_ptr);

    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &emissivity = p.get<MAST::FieldFunction<Real> >("emissivity");
    
    const MAST::Parameter
    &T_amb      = p.get<MAST::Parameter>("ambient_temperature"),
    &T_ref_zero = p.get<MAST::Parameter>("reference_zero_temperature"),
    &sb_const   = p.get<MAST::Parameter>("stefan_bolzmann_constant");
    
    
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
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location (qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        emissivity(pt, _time, emiss);
        temp  = phi_vec.dot(_sol);
        
        f   += JxW[qp] * phi_vec * sbc * emiss *
        (pow(temp-zero_ref, 4.) - pow(amb_temp-zero_ref, 4.));
        
        if (request_jacobian) {
            
            Bmat.reinit(1, phi_vec);
            Bmat.right_multiply_transpose(mat, Bmat);
            jac +=  JxW[qp] * mat * sbc * emiss * 4. * pow(temp-zero_ref, 3.);
        }
    }
    
    
    return request_jacobian;
}




bool
MAST::HeatConductionElementBase::
surface_radiation_residual(bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac,
                           MAST::BoundaryConditionBase& p) {
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>
    &emissivity = p.get<MAST::FieldFunction<Real> >("emissivity");
    
    const MAST::Parameter
    &T_amb      = p.get<MAST::Parameter>("ambient_temperature"),
    &T_ref_zero = p.get<MAST::Parameter>("reference_zero_temperature"),
    &sb_const   = p.get<MAST::Parameter>("stefan_bolzmann_constant");
    
    
    const std::vector<Real> &JxW               = _fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint  = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = _fe->get_phi();
    const unsigned int n_phi                   = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    RealMatrixX mat      = RealMatrixX::Zero(n_phi, n_phi);
    const Real
    sbc      = sb_const(),
    amb_temp = T_amb(),
    zero_ref = T_ref_zero();
    Real temp, emiss;
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location (qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // value of flux
        emissivity(pt, _time, emiss);
        temp  = phi_vec.dot(_sol);
        
        f   += JxW[qp] * phi_vec * sbc * emiss *
        (pow(temp-zero_ref, 4.) - pow(amb_temp-zero_ref, 4.));
        
        if (request_jacobian) {
            
            Bmat.reinit(1, phi_vec);
            Bmat.right_multiply_transpose(mat, Bmat);
            jac +=  JxW[qp] * mat * sbc * emiss * 4. * pow(temp-zero_ref, 3.);
        }
    }
    
    
    return request_jacobian;
}





bool
MAST::HeatConductionElementBase::
surface_radiation_residual_sensitivity(bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac,
                                       const unsigned int s,
                                       MAST::BoundaryConditionBase& p) {
    
    return true;
}




bool
MAST::HeatConductionElementBase::
surface_radiation_residual_sensitivity(bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac,
                                       MAST::BoundaryConditionBase& p) {
    
    return true;
}





bool
MAST::HeatConductionElementBase::
volume_heat_source_residual(bool request_jacobian,
                            RealVectorX& f,
                            RealMatrixX& jac,
                            MAST::BoundaryConditionBase& p) {
    
    
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>& func =
    p.get<MAST::FieldFunction<Real> >("heat_source");
    
    
    const std::vector<Real> &JxW                 = _fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint    = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi   = _fe->get_phi();
    const unsigned int n_phi                     = (unsigned int)phi.size();
    
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    libMesh::Point pt;
    Real  source;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        _local_elem->global_coordinates_location (qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        // get the value of heat source
        func(pt, _time, source);
        
        f   -=  JxW[qp] * phi_vec * source;
    }
    
    // calculation of the load vector is independent of solution
    return false;
}



bool
MAST::HeatConductionElementBase::
volume_heat_source_residual_sensitivity(bool request_jacobian,
                                        RealVectorX& f,
                                        RealMatrixX& jac,
                                        MAST::BoundaryConditionBase& p) {
    
    return false;
}



void
MAST::HeatConductionElementBase::
_initialize_mass_fem_operator(const unsigned int qp,
                              const libMesh::FEBase& fe,
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
                                  const libMesh::FEBase& fe,
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


