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
#include "elasticity/piston_theory_boundary_condition.h"



MAST::PistonTheoryPressure::
PistonTheoryPressure(unsigned int order,
                     MAST::FieldFunction<Real> *V,
                     MAST::FieldFunction<Real> *M,
                     MAST::FieldFunction<Real> *rho,
                     MAST::FieldFunction<Real> *gamma,
                     MAST::FieldFunction<Real> *dwdx,
                     MAST::FieldFunction<Real> *dwdt):
MAST::FieldFunction<Real>("pressure"),
_order(order),
_V_inf(V),
_M_inf(M),
_rho_inf(rho),
_gamma(gamma),
_dwdx(dwdx),
_dwdt(dwdt) {
    
    _functions.insert(_V_inf->master());
    _functions.insert(_M_inf->master());
    _functions.insert(_rho_inf->master());
    _functions.insert(_gamma->master());
    _functions.insert(_dwdx->master());
    _functions.insert(_dwdt->master());
}



MAST::PistonTheoryPressure::
PistonTheoryPressure(const MAST::PistonTheoryPressure &f):
MAST::FieldFunction<Real>(f),
_order(f._order),
_V_inf(f._V_inf->clone().release()),
_M_inf(f._M_inf->clone().release()),
_rho_inf(f._rho_inf->clone().release()),
_gamma(f._gamma->clone().release()),
_dwdx(f._dwdx->clone().release()),
_dwdt(f._dwdt->clone().release()) {
    
    libmesh_assert_less_equal(_order, 3);
    _functions.insert(_V_inf->master());
    _functions.insert(_M_inf->master());
    _functions.insert(_rho_inf->master());
    _functions.insert(_gamma->master());
    _functions.insert(_dwdx->master());
    _functions.insert(_dwdt->master());
}




std::auto_ptr<MAST::FieldFunction<Real> >
MAST::PistonTheoryPressure::clone() const {
    
    return std::auto_ptr<MAST::FieldFunction<Real> >
    (new MAST::PistonTheoryPressure(*this));
}



MAST::PistonTheoryPressure::~PistonTheoryPressure() {
    
    delete _V_inf;
    delete _M_inf;
    delete _rho_inf;
    delete _gamma;
}



void
MAST::PistonTheoryPressure::operator() (const libMesh::Point& p,
                                        const Real t,
                                        Real& m) const {
    
    Real V, M, rho, gamma, dwdx, dwdt;
    (*_V_inf)   (p, t,     V);
    (*_M_inf)   (p, t,     M);
    (*_rho_inf) (p, t,   rho);
    (*_gamma)   (p, t, gamma);
    (*_dwdx)    (p, t,  dwdx);
    (*_dwdt)    (p, t,  dwdt);
    
    
    Real
    v0  = dwdx + (M*M-2)/(M*M-1)*dwdt/V,
    v1  = V*V*dwdx + (M*M-2)/(M*M-1)*V*dwdt,
    v2  = pow(V*dwdx + (M*M-2)/(M*M-1)*dwdt, 2),
    v3  = v0*v2;
    
    m = v1; // linear
    switch (_order) {
            
        case 3:
            m += (gamma+1)/12.*M*M*v3;
        case 2:
            m += (gamma+1)/4.*M*v2;
    }
    
    m *= rho/sqrt(M*M-1.);
}




void
MAST::PistonTheoryPressure::derivative (const MAST::DerivativeType d,
                                        const MAST::FunctionBase& f,
                                        const libMesh::Point& p,
                                        const Real t,
                                        Real& m) const {
    
    Real V, M, rho, gamma, dwdx, dwdt, dV, dM, drho, dgamma, ddwdx, ddwdt;
    (*_V_inf)   (p, t,     V);  _V_inf->derivative   (d, f, p, t,     dV);
    (*_M_inf)   (p, t,     M);  _M_inf->derivative   (d, f, p, t,     dM);
    (*_rho_inf) (p, t,   rho);  _rho_inf->derivative (d, f, p, t,   drho);
    (*_gamma)   (p, t, gamma);  _gamma->derivative   (d, f, p, t, dgamma);
    (*_dwdx)    (p, t,  dwdx);  _dwdx->derivative    (d, f, p, t, ddwdx);
    (*_dwdt)    (p, t,  dwdt);  _dwdt->derivative    (d, f, p, t, ddwdt);
    
    Real
    Mf      = (M*M-2)/(M*M-1),
    dMf     = 2*M*dM/(M*M-1) - (M*M-2)/pow(M*M-1,2)*2*M*dM,
    v0      = dwdx + Mf*dwdt/V,
    v1      = V*V*dwdx + Mf*V*dwdt,
    v2      = pow(V*dwdx + Mf*dwdt, 2),
    v3      = v0*v2,
    v0dw    = ddwdx,
    v1dw    = V*V*ddwdx,
    v2dw    = 2.*(V*dwdx + Mf*dwdt)*V*ddwdx,
    v3dw    = v0*v2dw+v0dw*v2,
    v0dwt   = Mf/V*ddwdt,
    v1dwt   = Mf*V*ddwdt,
    v2dwt   = 2.*(V*dwdx + Mf*dwdt)*Mf*ddwdt,
    v3dwt   = v0dwt*v2+v0*v2dwt,
    dv0dV   = -Mf*dwdt/(V*V)*dV,
    dv1dV   = (2.*V*dwdx + Mf*dwdt)*dV,
    dv2dV   = 2.*(V*dwdx + Mf*dwdt)*dwdx*dV,
    dv3dV   = dv0dV*v2+v0*dv2dV,
    dv0dM   = dMf*dwdt/V,
    dv1dM   = dMf*V*dwdt,
    dv2dM   = 2.*(V*dwdx + Mf*dwdt)*dMf*dwdt,
    dv3dM   = dv0dM*v2+v0*dv2dM,
    f1      = rho/sqrt(M*M-1.),
    df1drho = drho/sqrt(M*M-1.),
    df1dM   = -.5*rho/pow(M*M-1., 1.5)*2*M*dM;
    
    
    //
    //  linear part:  f1 * v1,
    //  where    f1 = rho/sqrt(M*M-1.) and
    //           v1 = V*V*dwdx + Mf*V*dwdt,
    //
    m = f1*(dv1dV + dv1dM + v1dw + v1dwt) + (df1drho + df1dM) * v1; // linear
    
    switch (_order) {
            
        case 3:
            //
            //  cubic part:  f1 * (gamma+1)/12.*M*M*v3
            //
            m +=
            (df1drho + df1dM) *M*M*v3 * (gamma+1)/12.+
            f1*M*(2*dM*v3 + M*(dv3dV + dv3dM + v3dw + v3dwt)) * (gamma+1)/12. +
            f1 * M*M*v3/12.*dgamma;
        case 2:
            //
            //  quadratic part:  f1 * (gamma+1)/4.*M*v2
            //
            m +=
            (df1drho + df1dM) * M*v2 * (gamma+1)/4.+
            f1*(dM*v2+M*(dv2dV + dv2dM + v2dw + v2dwt)) * (gamma+1)/4. +
            f1*M*v2/4.*dgamma;
    }
}



MAST::PistonTheoryBoundaryCondition::
PistonTheoryBoundaryCondition(unsigned int order,
                              const RealVectorX& vel_vec):
MAST::BoundaryConditionBase(MAST::PISTON_THEORY),
_order(order),
_vel_vec(vel_vec) {
    
    // add the parameters
    
    
    // scale the velocity vector to unit length
    _vel_vec.normalize();
}


MAST::PistonTheoryBoundaryCondition::~PistonTheoryBoundaryCondition() {
    
}



unsigned int
MAST::PistonTheoryBoundaryCondition::order() const {
    return _order;
}



const RealVectorX&
MAST::PistonTheoryBoundaryCondition::vel_vec() const {
    return _vel_vec;
}



std::auto_ptr<MAST::PistonTheoryPressure>
MAST::PistonTheoryBoundaryCondition::
get_pressure_function(MAST::FieldFunction<Real>& dwdx,
                      MAST::FieldFunction<Real>& dwdt) const {
    
    MAST::PistonTheoryPressure
    *rval = new MAST::PistonTheoryPressure
    (_order,
     this->get<MAST::FieldFunction<Real> >("V").clone().release(),
     this->get<MAST::FieldFunction<Real> >("mach").clone().release(),
     this->get<MAST::FieldFunction<Real> >("rho").clone().release(),
     this->get<MAST::FieldFunction<Real> >("gamma").clone().release(),
     dwdx.clone().release(),
     dwdt.clone().release());
    
    return std::auto_ptr<MAST::PistonTheoryPressure>(rval);
}



