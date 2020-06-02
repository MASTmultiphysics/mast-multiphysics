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
#include "level_set/level_set_elem_base.h"
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
#include "property_cards/element_property_card_1D.h"


MAST::LevelSetElementBase::
LevelSetElementBase(MAST::SystemInitialization&             sys,
                    const MAST::GeomElem&                   elem):
MAST::ElementBase(sys, elem),
_phi_vel          (nullptr),
_if_propagation   (true) {
    
}



MAST::LevelSetElementBase::~LevelSetElementBase() {
    
}



void
MAST::LevelSetElementBase::
set_reference_solution_for_initialization(const RealVectorX& sol) {

    libmesh_assert(!_if_propagation);
    _ref_sol = sol;
}



bool
MAST::LevelSetElementBase::internal_residual (bool request_jacobian,
                                              RealVectorX& f,
                                              RealMatrixX& jac) {
    
    libmesh_assert(_phi_vel);
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(true, false,
                                                   _system.system().extra_quadrature_order));
    
    const std::vector<Real>& JxW           = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = fe->get_xyz();
    const unsigned int
    n_phi  = fe->n_shape_functions(),
    dim    = _elem.dim();
    
    RealMatrixX
    eye            = RealMatrixX::Identity(1, 1),
    tau            = RealMatrixX::Zero(    1, 1),
    mat1_n1n2      = RealMatrixX::Zero(    1, n_phi),
    mat2_n1n2      = RealMatrixX::Zero(    1, n_phi),
    mat_n2n2       = RealMatrixX::Zero(n_phi, n_phi);
    RealVectorX
    vec1_n1  = RealVectorX::Zero(1),
    vec2_n1  = RealVectorX::Zero(1),
    vec2_n2  = RealVectorX::Zero(n_phi),
    flux     = RealVectorX::Zero(1),
    vel      = RealVectorX::Zero(dim);
    Real
    dc       = 0.,
    source   = 0.;
    
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix Bmat;
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_fem_operators(qp, *fe, Bmat, dBmat);
        
        _velocity_and_source(qp, xyz[qp], _time, Bmat, dBmat, vel, source);
        _tau(*fe, qp, Bmat, dBmat, vel, tau);
        _dc_operator(*fe, qp, dBmat, vel, dc);
        
        // calculate the flux for each dimension and add its weighted
        // component to the residual
        flux.setZero();
        mat2_n1n2.setZero();
        for (unsigned int j=0; j<dim; j++) {
            
            dBmat[j].right_multiply(vec1_n1, _sol);       // dphi_dx_j
            flux += vel(j) * vec1_n1;                     // v_j dphi_dx_j
            dBmat[j].left_multiply(mat1_n1n2, eye);       // dBmat
            mat2_n1n2 += mat1_n1n2 * vel(j);              // dBmat V_j
        }

        // add to the residual vector
        flux(0) -= source;                                   // V.grad(phi) - s
        Bmat.vector_mult_transpose(vec2_n2, flux);
        f += JxW[qp] * vec2_n2;                              // int_omega          u      (V.grad(phi)-s)
        f += JxW[qp] * mat2_n1n2.transpose() * tau * flux;   // int_omega   V.grad(u) tau (V.grad(phi)-s)
        
        // discontinuity capturing
        for (unsigned int j=0; j<dim; j++) {
            
            dBmat[j].vector_mult(vec1_n1, _sol);
            dBmat[j].vector_mult_transpose(vec2_n2, vec1_n1);
            f += JxW[qp] * dc * vec2_n2;
        }

        
        if (request_jacobian) {
            
            for (unsigned int j=0; j<dim; j++) {
                    
                Bmat.right_multiply_transpose(mat_n2n2, dBmat[j]);
                jac += JxW[qp] * vel(j) * mat_n2n2;  // int_omega  u  V.grad(phi)
                
                // discontinuity capturing term
                dBmat[j].right_multiply_transpose(mat_n2n2, dBmat[j]);   // dB_i^T dc dB_i
                jac += JxW[qp] * dc * mat_n2n2;
            }
            
            jac += JxW[qp] * mat2_n1n2.transpose() * tau * mat2_n1n2;  // int_omega  V.grad(u) tau (V.grad(phi))
        }
    }
    
    return request_jacobian;
}





bool
MAST::LevelSetElementBase::velocity_residual (bool request_jacobian,
                                              RealVectorX& f,
                                              RealMatrixX& jac_xdot,
                                              RealMatrixX& jac) {

    libmesh_assert(_phi_vel);

    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(true, false,
                                                   _system.system().extra_quadrature_order));

    const std::vector<Real>& JxW           = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = fe->get_xyz();
    const unsigned int
    n_phi  = fe->n_shape_functions(),
    dim    = _elem.dim();
    
    RealMatrixX
    eye            = RealMatrixX::Identity(1, 1),
    tau            = RealMatrixX::Zero(    1, 1),
    mat1_n1n2      = RealMatrixX::Zero(    1, n_phi),
    mat2_n1n2      = RealMatrixX::Zero(    1, n_phi),
    mat_n2n1       = RealMatrixX::Zero(n_phi,     1),
    mat_n2n2       = RealMatrixX::Zero(n_phi, n_phi);
    RealVectorX
    vec1_n1  = RealVectorX::Zero(1),
    vec2_n2  = RealVectorX::Zero(n_phi),
    flux     = RealVectorX::Zero(1),
    vel      = RealVectorX::Zero(dim);
    Real
    source   = 0.;

    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix Bmat;
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_fem_operators(qp, *fe, Bmat, dBmat);
        
        _velocity_and_source(qp, xyz[qp], _time, Bmat, dBmat, vel, source);
        _tau(*fe, qp, Bmat, dBmat, vel, tau);
        
        // calculate the flux for each dimension and add its weighted
        // component to the residual
        mat2_n1n2.setZero();
        for (unsigned int j=0; j<dim; j++) {
            
            dBmat[j].left_multiply(mat1_n1n2, eye);       // dBmat
            mat2_n1n2 += mat1_n1n2 * vel(j);              // dBmat_j V_j
        }
        
        // add to the residual vector
        Bmat.right_multiply(vec1_n1, _vel);                     // dphi/dt
        Bmat.vector_mult_transpose(vec2_n2, vec1_n1);
        f += JxW[qp] * vec2_n2;                                 // int_omega          u      dphi/dt
        f += JxW[qp] * mat2_n1n2.transpose() * tau * vec1_n1;   // int_omega   V.grad(u) tau dphi/dt
        
        if (request_jacobian) {
            
            Bmat.right_multiply_transpose(mat_n2n2, Bmat);
            jac_xdot += JxW[qp] * mat_n2n2;  // int_omega         u      dphi/dt
            
            mat_n2n1  = mat2_n1n2.transpose()*tau;
            Bmat.left_multiply(mat_n2n2, mat_n2n1);
            jac_xdot += JxW[qp] * mat_n2n2;  // int_omega  V.grad(u) tau dphi/dt
        }
    }
    
    return request_jacobian;
}




bool
MAST::LevelSetElementBase::
side_external_residual (bool request_jacobian,
                        RealVectorX& f,
                        RealMatrixX& jac,
                        std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    libmesh_assert(false);
    return false;
}





bool
MAST::LevelSetElementBase::
volume_external_residual (bool request_jacobian,
                          RealVectorX& f,
                          RealMatrixX& jac,
                          std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    libmesh_assert(false);
    return false;
}



bool
MAST::LevelSetElementBase::
side_external_residual_sensitivity (bool request_jacobian,
                                    RealVectorX& f,
                                    RealMatrixX& jac,
                                    std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    libmesh_assert(false);
    return false;
}




bool
MAST::LevelSetElementBase::
volume_external_residual_sensitivity (bool request_jacobian,
                                      RealVectorX& f,
                                      RealMatrixX& jac,
                                      std::multimap<libMesh::subdomain_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    libmesh_assert(false);
    return false;
}




bool
MAST::LevelSetElementBase::
internal_residual_sensitivity (bool request_jacobian,
                               RealVectorX& f,
                               RealMatrixX& jac) {
    
    libmesh_assert(false); // should not get called
    return request_jacobian;
}



bool
MAST::LevelSetElementBase::
velocity_residual_sensitivity (bool request_jacobian,
                               RealVectorX& f,
                               RealMatrixX& jac) {
    
    libmesh_assert(false); // should not get called.
    return request_jacobian;
}



Real
MAST::LevelSetElementBase::volume() {
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(true, false,
                                                   _system.system().extra_quadrature_order));

    const std::vector<Real>& JxW           = fe->get_JxW();
    const unsigned int
    dim    = _elem.dim();
    
    RealVectorX
    phi      = RealVectorX::Zero(1);

    Real
    vol      = 0.;
    
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix Bmat;
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_fem_operators(qp, *fe, Bmat, dBmat);
        Bmat.right_multiply(phi, _sol);
        if (phi(0) > 0.) vol += JxW[qp];
    }
    
    return vol;
}


Real
MAST::LevelSetElementBase::homogenized_volume_fraction(Real delta) {
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(true, false,
                                                   _system.system().extra_quadrature_order));
    
    const std::vector<Real>& JxW           = fe->get_JxW();
    const unsigned int
    dim    = _elem.dim();
    
    RealVectorX
    phi      = RealVectorX::Zero(1);
    
    Real
    pi       = acos(-1.),
    heav     = 0.,
    vol      = 0.;
    
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix Bmat;
    

    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_fem_operators(qp, *fe, Bmat, dBmat);
        Bmat.right_multiply(phi, _sol);
        vol  += JxW[qp];
        heav += .5 * (1. + 2./pi * atan(phi(0)/delta)) * JxW[qp];
    }
    
    return heav/vol;
}



Real
MAST::LevelSetElementBase::homogenized_volume_fraction_sensitivity(Real delta) {
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(true, false,
                                                   _system.system().extra_quadrature_order));
    
    const std::vector<Real>& JxW           = fe->get_JxW();
    const unsigned int
    dim    = _elem.dim();
    
    RealVectorX
    phi      = RealVectorX::Zero(1),
    dphidp   = RealVectorX::Zero(1);
    
    Real
    pi       = acos(-1.),
    heav     = 0.,
    vol      = 0.;
    
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix Bmat;
    

    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_fem_operators(qp, *fe, Bmat, dBmat);
        Bmat.right_multiply(phi,         _sol);
        Bmat.right_multiply(dphidp, _sol_sens);
        vol  += JxW[qp];
        heav += 1./pi/delta/(1.+pow(phi(0)/delta, 2)) * dphidp(0) * JxW[qp];
    }
    
    return heav/vol;
}



Real
MAST::LevelSetElementBase::perimeter(Real d) {
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(true, false,
                                                   _system.system().extra_quadrature_order));
    
    const std::vector<Real>& JxW           = fe->get_JxW();
    const unsigned int
    dim    = _elem.dim();
    
    RealVectorX
    phi      = RealVectorX::Zero(1);
    
    Real
    pi       = acos(-1.),
    per      = 0.;
    
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix Bmat;
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_fem_operators(qp, *fe, Bmat, dBmat);
        Bmat.right_multiply(phi, _sol);
        per += 1./pi/d/(1.+pow(phi(0)/d, 2)) * JxW[qp];
    }
    
    return per;
}


Real
MAST::LevelSetElementBase::perimeter_sensitivity(Real d) {
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_fe(true, false,
                                                   _system.system().extra_quadrature_order));
    
    const std::vector<Real>& JxW           = fe->get_JxW();
    const unsigned int
    dim    = _elem.dim();
    
    RealVectorX
    phi      = RealVectorX::Zero(1),
    dphidp   = RealVectorX::Zero(1);
    
    Real
    pi       = acos(-1.),
    dper_dp  = 0.;
    
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix Bmat;
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_fem_operators(qp, *fe, Bmat, dBmat);
        Bmat.right_multiply(phi,         _sol);
        Bmat.right_multiply(dphidp, _sol_sens);
        dper_dp -= 2.*phi(0)/pi/pow(d,3)/pow(1.+pow(phi(0)/d, 2),2) * dphidp(0) * JxW[qp];
    }
    
    return dper_dp;
}




Real
MAST::LevelSetElementBase::volume_boundary_velocity_on_side(unsigned int s) {
    
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(s, true, false));

    const std::vector<Real>& JxW           =  fe->get_JxW();
    const unsigned int
    dim    = _elem.dim();
    
    RealVectorX
    phi      = RealVectorX::Zero(1),
    gradphi  = RealVectorX::Zero(dim),
    dphidp   = RealVectorX::Zero(1),
    vel      = RealVectorX::Zero(dim);
    
    Real
    vn       = 0.,
    dvoldp   = 0.;
    
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix Bmat;
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_fem_operators(qp, *fe, Bmat, dBmat);

        // first calculate gradient of phi
        for (unsigned int i=0; i<dim; i++) {
            dBmat[i].right_multiply(phi, _sol);
            gradphi(i) = phi(0);
        }

        // initialize the Bmat operator for this term
        _initialize_fem_operators(qp, *fe, Bmat, dBmat);
        Bmat.right_multiply(phi,         _sol);
        Bmat.right_multiply(dphidp, _sol_sens);
        
        // at boundary, phi(x) = 0
        // so,  dphi/dp + grad(phi) . V = 0
        //      dphi/dp + grad(phi) . (-grad(phi)/|grad(phi)|  Vn) = 0   [since V = -grad(phi)/|grad(phi)| Vn]
        //      dphi/dp -(grad(phi) .  grad(phi))/|grad(phi)|  Vn  = 0
        //      Vn  =  (dphi/dp) |grad(phi)|  / |grad(phi)|^2 = (dphi/dp) / |grad(phi)|
        vn      =  dphidp(0) / gradphi.norm();

        // vol     =  int_omega  dOmega
        // dvol/dp =  int_Gamma Vn dGamma
        dvoldp += vn * JxW[qp];
    }
    
    return dvoldp;
}



void
MAST::LevelSetElementBase::_velocity_and_source(const unsigned int qp,
                                                const libMesh::Point& p,
                                                const Real t,
                                                const MAST::FEMOperatorMatrix& Bmat,
                                                const std::vector<MAST::FEMOperatorMatrix>& dBmat,
                                                RealVectorX& vel,
                                                Real&        source) {
    
    const unsigned int
    dim       =  _elem.dim();
    
    RealVectorX
    v         =  RealVectorX::Zero(1),
    grad_phi  =  RealVectorX::Zero(dim);
    
    // first initialize grad(phi)
    for (unsigned int i=0; i<dim; i++) {
        
        dBmat[i].right_multiply(v, _sol);
        grad_phi(i) = v(0);
    }
    
    if (_if_propagation) {

        // now, initialize the velocity vector
        Real
        Vn        = 0.;
        
        (*_phi_vel)(p, t, Vn);
        vel       =  grad_phi * (-Vn/grad_phi.norm());
        source    = 0.;
    }
    else {
        
        libmesh_assert_equal_to(_ref_sol.size(), _sol.size());
        
        Bmat.right_multiply(v, _ref_sol);
        
        Real
        tol           = 1.e-6,
        ref_phi       = v(0),
        grad_phi_norm = grad_phi.norm();
        
        source = ref_phi/pow(pow(ref_phi,2)+ tol*pow(grad_phi_norm,2), 0.5);
        
        if (grad_phi_norm > tol)
            vel        = grad_phi * source/grad_phi.norm();
        else
            vel.setZero();
    }
}


void
MAST::LevelSetElementBase::_tau(const MAST::FEBase& fe,
                                unsigned int qp,
                                const MAST::FEMOperatorMatrix& Bmat,
                                const std::vector<MAST::FEMOperatorMatrix>& dBmat,
                                const RealVectorX& vel,
                                RealMatrixX& tau) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe.get_dphi();
    const unsigned int n_phi = (unsigned int)fe.n_shape_functions();
    RealVectorX
    phi = RealVectorX::Zero(n_phi);

    Real
    tol    = 1.e-6,
    val1   = 0.,
    val2   = 0.,
    vel_l2 = vel.norm();

    libMesh::Point
    nvec;
    
    if (vel_l2 > tol) {
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ ) {
         
            nvec  = dphi[i_nd][qp];
            val2 = 0.;
            for (unsigned int i_dim=0; i_dim<_elem.dim(); i_dim++)
                val2 += nvec(i_dim) * vel(i_dim);
            
            val1 += std::fabs(val2);
        }
        
        tau(0,0) = 1./val1;
    }
    else
        tau(0,0) = 0.;
}




void
MAST::LevelSetElementBase::_dc_operator(const MAST::FEBase& fe,
                                        const unsigned int qp,
                                        const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
                                        const RealVectorX& vel,
                                        Real& dc) {

    unsigned int
    dim = _elem.dim();

    RealVectorX
    dphi               = RealVectorX::Zero(dim),
    vec1               = RealVectorX::Zero(1),
    dflux              = RealVectorX::Zero(1);
    
    RealMatrixX
    dxi_dX             = RealMatrixX::Zero(dim, dim);

    
    _calculate_dxidX(fe, qp, dxi_dX);

    for (unsigned int i=0; i<dim; i++) {

        dB_mat[i].vector_mult(vec1, _sol); // dphi/dx_i
        dflux(0) += vel(i) * vec1(0);      // V.grad(phi)
        dphi(i)   = vec1(0);               // grad(phi)
    }

    Real
    val1 = 1.0e-6;
    
    for (unsigned int i=0; i<dim; i++) {

        vec1.setZero();
        
        for (unsigned int j=0; j<dim; j++)
            vec1(0) += dxi_dX(i, j) * dphi(j);
        
        // calculate the value of denominator
        val1 += pow(vec1.norm(), 2);
    }

    dc = 0.5*dflux.norm()/pow(val1, 0.5);
}



void
MAST::LevelSetElementBase::
_calculate_dxidX (const MAST::FEBase& fe,
                  const unsigned int qp,
                  RealMatrixX& dxi_dX) {
    
    
    // initialize dxi_dX and dX_dxi
    unsigned int
    dim   = _elem.dim();
    Real
    val   = 0.;
    dxi_dX.setZero();

    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        for (unsigned int j_dim=0; j_dim<dim; j_dim++)
        {
            switch (i_dim)
            {
                case 0:
                {
                    switch (j_dim)
                    {
                        case 0:
                            val  = fe.get_dxidx()[qp];
                            break;
                        case 1:
                            val = fe.get_dxidy()[qp];
                            break;
                        case 2:
                            val = fe.get_dxidz()[qp];
                            break;
                    }
                }
                    break;
                case 1:
                {
                    switch (j_dim)
                    {
                        case 0:
                            val = fe.get_detadx()[qp];
                            break;
                        case 1:
                            val = fe.get_detady()[qp];
                            break;
                        case 2:
                            val = fe.get_detadz()[qp];
                            break;
                    }
                }
                    break;
                case 2:
                {
                    switch (j_dim)
                    {
                        case 0:
                            val = fe.get_dzetadx()[qp];
                            break;
                        case 1:
                            val = fe.get_dzetady()[qp];
                            break;
                        case 2:
                            val = fe.get_dzetadz()[qp];
                            break;
                    }
                }
                    break;
            }
            dxi_dX(i_dim, j_dim) = val;
        }
}



void
MAST::LevelSetElementBase::
_initialize_fem_operators(const unsigned int qp,
                          const MAST::FEBase& fe,
                          MAST::FEMOperatorMatrix& Bmat,
                          std::vector<MAST::FEMOperatorMatrix>& dBmat) {
    
    const std::vector<std::vector<Real> >&                   phi_fe = fe.get_phi();
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe.get_dphi();

    const unsigned int n_phi = (unsigned int)phi_fe.size();
    
    RealVectorX
    phi = RealVectorX::Zero(n_phi);
    
    // shape function values
    // N
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = phi_fe[i_nd][qp];
    
    Bmat.reinit(1, phi);
    
    // now set the shape function derivatives
    for (unsigned int i_dim=0; i_dim<_elem.dim(); i_dim++) {

        phi.setZero();
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi(i_nd) = dphi[i_nd][qp](i_dim);
        dBmat[i_dim].reinit(1, phi); //  dT/dx_i
    }

}


