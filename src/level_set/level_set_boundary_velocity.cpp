/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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
#include "level_set/level_set_boundary_velocity.h"


MAST::LevelSetBoundaryVelocity::LevelSetBoundaryVelocity(const unsigned int dim):
MAST::FieldFunction<RealVectorX>("phi_vel"),
_dim(dim),
_phi(nullptr) {
    
}

MAST::LevelSetBoundaryVelocity::~LevelSetBoundaryVelocity() {
    
    if (_phi)
        delete _phi;
}


void
MAST::LevelSetBoundaryVelocity::init(MAST::SystemInitialization& sys,
                                     const libMesh::NumericVector<Real>& sol,
                                     const libMesh::NumericVector<Real>* dsol) {
    
    if (!_phi)
        _phi = new MAST::MeshFieldFunction(sys, "phi_vel");
    else
        _phi->clear();
    _phi->init(sol, dsol);
}


void
MAST::LevelSetBoundaryVelocity::operator() (const libMesh::Point& p,
                                            const Real t,
                                            RealVectorX& v) const {
    this->velocity(p, t, v);
}


void
MAST::LevelSetBoundaryVelocity::velocity(const libMesh::Point& p,
                                         const Real t,
                                         RealVectorX& v) const {
    
    libmesh_assert(_phi);
    
    RealVectorX
    val      = RealVectorX::Zero(1),
    grad     = RealVectorX::Zero(_dim),
    dval     = RealVectorX::Zero(1);
    
    RealMatrixX
    gradmat  = RealMatrixX::Zero(1, _dim);
    
    // the velocity is identified using the level set function gradient
    // and its sensitivity
    (*_phi)(p, t, val);
    
    // since the function provides the velocity at the boundary, the
    // level set value should be close to zero. This is commented out since
    // for a coarse level set mesh used to define geometry on a fine analysis
    // mesh, the boundary can exist on locations that are not necessarily phi=0.
    //libmesh_assert_less_equal(std::fabs(val(0)), tol);
    //if (std::fabs(val(0)) > tol)
    //    libMesh::out << "**** !!!  level-set not zero at boundary: " << val(0) << std::endl;
    
    _phi->gradient(p, t, gradmat);
    _phi->perturbation(p, t, dval);
    
    // copy the matrix to the vector for further calculations
    grad = gradmat.row(0);
    
    // at boundary, phi(x) = 0
    // so,  dphi/dp + grad(phi) . V = 0
    //      dphi/dp + grad(phi) . (-grad(phi)/|grad(phi)|  Vn) = 0   [since V = -grad(phi)/|grad(phi)| Vn]
    //      dphi/dp -(grad(phi) .  grad(phi))/|grad(phi)|  Vn  = 0
    //      Vn  =  (dphi/dp)   |grad(phi)|  / |grad(phi)|^2
    //          =  (dphi/dp) / |grad(phi)|
    //      V   =  -grad(phi)/ |grad(phi)| Vn
    //          =  -grad(phi)/ |grad(phi)| * (dphi/dp) / |grad(phi)|
    //          =  -grad(phi) * (dphi/dp)/|grad(phi)|^2
    //
    v = -dval(0)/grad.squaredNorm()*grad;
}




void
MAST::LevelSetBoundaryVelocity::
search_nearest_interface_point(const libMesh::Point& p,
                               const Real t,
                               const Real length,
                               RealVectorX& v) const {
    
    libmesh_assert(_phi);
    
    RealVectorX
    phi      = RealVectorX::Zero(1),
    grad_phi = RealVectorX::Zero(_dim),
    dval     = RealVectorX::Zero(1),
    p_ref    = RealVectorX::Zero(3),
    dv0      = RealVectorX::Zero(4),
    dv1      = RealVectorX::Zero(4),
    gradL    = RealVectorX::Zero(4);
    
    RealMatrixX
    gradmat  = RealMatrixX::Zero(1, _dim),
    coeffs   = RealMatrixX::Zero(3, 3);

    libMesh::Point
    p_opt = p;
    
    p_ref(0) = p(0);
    p_ref(1) = p(1);
    p_ref(2) = p(2);

    // The point is identified from a constrained optimization problem.
    //      p*  = argmin { .5 * || p* - p ||^2:  phi(p*) = 0 }
    //  This is formulated as a constrained minimization problem:
    //      (p*, l)  = argmin { .5 * || p* - p ||^2 + l * phi(p*) },
    //  where the Lagrangian is defined as
    //      L (p*, l) = .5 * || p* - p ||^2 + l * phi(p*)
    //                = .5 * ((px* - px)^2 + (py* - py)^2 + (pz* - pz)^2) + l * phi(p*)
    //   The first order optimality for this provides
    //      grad ( L (p*, l)) = { (p*-p) + l * grad(phi(p*));  phi(p*) } = 0
    //   We solve this using Newton-Raphson with an approximate Hessian
    //
    
    bool
    if_cont = true;

    unsigned int
    n_iters   = 0,
    max_iters = 100;
    
    Real
    tol  = 1.e-8,
    L0   = 1.e12,
    damp = 0.6,
    L1   = 0.;
    
    // initialize the design point
    dv0.topRows(3) = p_ref;
    dv0(3) = 1.; // arbitrary value of Lagrange multiplier
    
    while (if_cont) {
        
        // compute the gradient of Lagrangian
        (*_phi)        (p_opt, t,     phi);
        _phi->gradient (p_opt, t, gradmat);

        // this defines the gradient of Lagrangian
        gradL.topRows(3) = (dv0.topRows(3)-p_ref) + dv0(3) * gradmat.row(0).transpose();
        gradL(3)         = phi(0);
        
        // approximate hessian
        coeffs.setZero(4,4);
        for (unsigned int i=0; i<3; i++) {
            coeffs(i,i) = 1.;
            coeffs(3,i) = coeffs(i,3) = gradmat(0,i);
        }

        Eigen::FullPivLU<RealMatrixX> solver(coeffs);
        dv0 -= damp*solver.solve(gradL);
        // update the design points and check for convergence
        p_opt(0) = dv0(0);
        p_opt(1) = dv0(1);
        p_opt(2) = dv0(2);
        
        L1        = _evaluate_point_search_obj(p, t, dv0);

        if (n_iters == max_iters) {
            
            libMesh::Point dp = p_opt - p;
            if_cont = false;
            libMesh::out
            << "Warning: nearest interface point search did not converge."
            << std::endl;
            libMesh::out
            << "  given pt: ("
            << p(0) << " , " << p(1) << " , " << p(2)
            << ") -> mapped pt: ("
            << p_opt(0) << " , " << p_opt(1) << " , " << p_opt(2)
            << ").  phi = " << phi(0)
            << "  d/h =  " << dp.norm()/length << std::endl;
        }
        if (std::fabs(L1) <= tol) if_cont = false;
        if (std::fabs((L1-L0)/L0) <= tol) if_cont = false;
        
        L0 = L1;
        n_iters++;
    }
    
    v.setZero();
    v.topRows(_dim) = dv0.topRows(_dim);
}



void
MAST::LevelSetBoundaryVelocity::
search_nearest_interface_point_derivative(const MAST::FunctionBase& f,
                                          const libMesh::Point& p,
                                          const Real t,
                                          const Real length,
                                          RealVectorX& v) const {
    
    libmesh_assert(_phi);

    // velocity at this interface point is given by the normal velocity
    // So, we first find the point and then compute its velocity
    
    this->search_nearest_interface_point(p, t, length, v);
    
    libMesh::Point
    pt(v(0), v(1), v(2));
    
    // now compute the velocity at this point.
    v.setZero();
    this->velocity(pt, t, v);
}




void
MAST::LevelSetBoundaryVelocity::normal_at_point(const libMesh::Point& p,
                                                const Real t,
                                                RealVectorX& n) const {
    
    libmesh_assert(_phi);
    
    RealMatrixX
    gradmat  = RealMatrixX::Zero(1, _dim);
    
    n.setZero();
    
    // the velocity is identified using the level set function gradient
    // and its sensitivity
    _phi->gradient(p, t, gradmat);
    n.topRows(3) = -gradmat.row(0);
    
    n /= n.norm();
}



void
MAST::LevelSetBoundaryVelocity::normal_derivative_at_point(const MAST::FunctionBase& f,
                                                           const libMesh::Point& p,
                                                           const Real t,
                                                           RealVectorX& n) const {
    
    libmesh_assert(_phi);
    
    RealMatrixX
    gradmat  = RealMatrixX::Zero(1, _dim),
    dgrad    = RealMatrixX::Zero(1, _dim);
    
    n.setZero();
    
    // the velocity is identified using the level set function gradient
    // and its sensitivity
    _phi->gradient(p, t, gradmat);
    _phi->perturbation_gradient(p, t, dgrad);
    
    RealVectorX
    v  = gradmat.row(0),
    dv = dgrad.row(0);
    
    n.topRows(3) = (-dv/v.norm() +
                    + 0.5/std::pow(v.norm(),1.5)* v.dot(dv) * v);
}



Real
MAST::LevelSetBoundaryVelocity::
_evaluate_point_search_obj(const libMesh::Point& p,
                           const Real t,
                           const RealVectorX& dv) const {
    
    libMesh::Point p1(dv(0), dv(1), dv(2));
    RealVectorX phi;
    (*_phi)(p1, t,     phi);
    p1 -= p;
    return .5 * p1.norm_sq() + dv(3)*phi(0);
}


