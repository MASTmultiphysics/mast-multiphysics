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
#include "elasticity/piston_theory_boundary_condition.h"



MAST::PistonTheoryPressure::
PistonTheoryPressure(unsigned int order,
                     const MAST::FieldFunction<Real>& V,
                     const MAST::FieldFunction<Real>& M,
                     const MAST::FieldFunction<Real>& rho,
                     const MAST::FieldFunction<Real>& gamma,
                     const MAST::FieldFunction<Real>& dwdx,
                     const MAST::FieldFunction<Real>& dwdt):
MAST::FieldFunction<Real>("pressure"),
_order(order),
_V_inf(V),
_M_inf(M),
_rho_inf(rho),
_gamma(gamma),
_dwdx(dwdx),
_dwdt(dwdt) {
    
    _functions.insert(&_V_inf);
    _functions.insert(&_M_inf);
    _functions.insert(&_rho_inf);
    _functions.insert(&_gamma);
    _functions.insert(&_dwdx);
    _functions.insert(&_dwdt);
}




MAST::PistonTheoryPressure::~PistonTheoryPressure() { }



void
MAST::PistonTheoryPressure::operator() (const libMesh::Point& p,
                                        const Real t,
                                        Real& m) const {
    
    Real V, M, rho, gamma, dwdx, dwdt;
    _V_inf   (p, t,     V);
    _M_inf   (p, t,     M);
    _rho_inf (p, t,   rho);
    _gamma   (p, t, gamma);
    _dwdx    (p, t,  dwdx);
    _dwdt    (p, t,  dwdt);
    
    
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
MAST::PistonTheoryPressure::derivative (       const MAST::FunctionBase& f,
                                        const libMesh::Point& p,
                                        const Real t,
                                        Real& m) const {
    
    Real V, M, rho, gamma, dwdx, dwdt, dV, dM, drho, dgamma, ddwdx, ddwdt;
    _V_inf   (p, t,     V);  _V_inf.derivative   (f, p, t,     dV);
    _M_inf   (p, t,     M);  _M_inf.derivative   (f, p, t,     dM);
    _rho_inf (p, t,   rho);  _rho_inf.derivative (f, p, t,   drho);
    _gamma   (p, t, gamma);  _gamma.derivative   (f, p, t, dgamma);
    _dwdx    (p, t,  dwdx);  _dwdx.derivative    (f, p, t, ddwdx);
    _dwdt    (p, t,  dwdt);  _dwdt.derivative    (f, p, t, ddwdt);
    
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





MAST::PistonTheoryPressureXDerivative::
PistonTheoryPressureXDerivative(unsigned int order,
                                const MAST::FieldFunction<Real>& V,
                                const MAST::FieldFunction<Real>& M,
                                const MAST::FieldFunction<Real>& rho,
                                const MAST::FieldFunction<Real>& gamma,
                                const MAST::FieldFunction<Real>& dwdx,
                                const MAST::FieldFunction<Real>& dwdt):
MAST::FieldFunction<Real>("pressure"),
_order(order),
_V_inf(V),
_M_inf(M),
_rho_inf(rho),
_gamma(gamma),
_dwdx(dwdx),
_dwdt(dwdt) {
    
    _functions.insert(&_V_inf);
    _functions.insert(&_M_inf);
    _functions.insert(&_rho_inf);
    _functions.insert(&_gamma);
    _functions.insert(&_dwdx);
    _functions.insert(&_dwdt);
}





MAST::PistonTheoryPressureXDerivative::~PistonTheoryPressureXDerivative() { }



void
MAST::PistonTheoryPressureXDerivative::operator() (const libMesh::Point& p,
                                                   const Real t,
                                                   Real& m) const {
    
    Real V, M, rho, gamma, dwdx, dwdt;
    _V_inf    (p, t,     V);
    _M_inf    (p, t,     M);
    _rho_inf  (p, t,   rho);
    _gamma    (p, t, gamma);
    _dwdx     (p, t,  dwdx);
    _dwdt     (p, t,  dwdt);
    
    
    Real
    Mf      = (M*M-2)/(M*M-1),
    
    v0      = dwdx + (M*M-2)/(M*M-1)*dwdt/V,
    v2      = pow(V*dwdx + (M*M-2)/(M*M-1)*dwdt, 2),
    
    v0dw    = 1.,
    v1dw    = V*V,
    v2dw    = 2.*(V*dwdx + Mf*dwdt)*V,
    v3dw    = v0*v2dw+v0dw*v2,
    
    f1      = rho/sqrt(M*M-1.);

    
    //
    //  linear part:  f1 * v1,
    //  where    f1 = rho/sqrt(M*M-1.) and
    //           v1 = V*V*dwdx + Mf*V*dwdt,
    //
    m = f1*v1dw; // linear
    
    switch (_order) {
            
        case 3:
            //
            //  cubic part:  f1 * (gamma+1)/12.*M*M*v3
            //
            m += f1*M*M*v3dw * (gamma+1)/12.;
        case 2:
            //
            //  quadratic part:  f1 * (gamma+1)/4.*M*v2
            //
            m += f1*M*v2dw * (gamma+1)/4.;
    }
}




void
MAST::PistonTheoryPressureXDerivative::derivative (                  const MAST::FunctionBase& f,
                                                   const libMesh::Point& p,
                                                   const Real t,
                                                   Real& m) const {
    
    Real V, M, rho, gamma, dwdx, dwdt, dV, dM, drho, dgamma;
    _V_inf    (p, t,     V);  _V_inf.derivative   (f, p, t,     dV);
    _M_inf    (p, t,     M);  _M_inf.derivative   (f, p, t,     dM);
    _rho_inf  (p, t,   rho);  _rho_inf.derivative (f, p, t,   drho);
    _gamma    (p, t, gamma);  _gamma.derivative   (f, p, t, dgamma);
    _dwdx     (p, t,  dwdx);
    _dwdt     (p, t,  dwdt);
    
    Real
    Mf      = (M*M-2)/(M*M-1),
    dMfdM   = 2*M*dM/(M*M-1) - (M*M-2)/pow(M*M-1,2)*2*M*dM,

    v0      = dwdx + Mf*dwdt/V,
    v2      = pow(V*dwdx + Mf*dwdt, 2),
    
    v0dw    = 1.,
    v1dw    = V*V,
    v2dw    = 2.*(V*dwdx + Mf*dwdt)*V,
    v3dw    = v0*v2dw+v0dw*v2,

    dv0dV   = -Mf*dwdt/(V*V)*dV,
    dv2dV   = 2.*(V*dwdx + Mf*dwdt)*dwdx*dV,
    
    dv0dM   = dMfdM*dwdt/V,
    dv2dM   = 2.*(V*dwdx + Mf*dwdt)*dMfdM*dwdt,

    v0dwdV  = 0.,
    v1dwdV  = 2.*V*dV,
    v2dwdV  = 2.*(2.*V*dV*dwdx + Mf*dV*dwdt),
    v3dwdV  = dv0dV*v2dw+v0*v2dwdV + v0dwdV*v2+v0dw*dv2dV,
    
    v0dwdM  = 0.,
    v1dwdM  = 0.,
    v2dwdM  = 2.*dMfdM*dwdt*V,
    v3dwdM  = dv0dM*v2dw+v0*v2dwdM + v0dwdM*v2+v0dw*dv2dM,
    
    f1      = rho/sqrt(M*M-1.),
    df1drho = drho/sqrt(M*M-1.),
    df1dM   = -.5*rho/pow(M*M-1., 1.5)*2*M*dM;
    
    
    //
    //  linear part:  d/dp (f1 * dv1/dw),
    //  where    f1 = rho/sqrt(M*M-1.) and
    //        dv1dw = V*V*dwdx,
    //
    m = f1*(v1dwdV + v1dwdM) + (df1drho + df1dM) * v1dw; // linear
    
    switch (_order) {
            
        case 3:
            //
            //  cubic part:  d/dp (f1 * (gamma+1)/12.*M*M* dv3/dw)
            //
            m +=
            (df1drho + df1dM) *M*M*v3dw * (gamma+1)/12.+
            f1*M* (2*dM*v3dw + M*(v3dwdV + v3dwdM)) * (gamma+1)/12. +
            f1 * M*M*v3dw/12.*dgamma;
        case 2:
            //
            //  quadratic part:  d/dp (f1 * (gamma+1)/4.*M* dv2/dw)
            //
            m +=
            (df1drho + df1dM) * M*v2dw * (gamma+1)/4.+
            f1*(dM*v2dw + M*(v2dwdV + v2dwdM)) * (gamma+1)/4. +
            f1*M*v2dw/4.*dgamma;
    }
}




MAST::PistonTheoryPressureXdotDerivative::
PistonTheoryPressureXdotDerivative(unsigned int order,
                                   const MAST::FieldFunction<Real>& V,
                                   const MAST::FieldFunction<Real>& M,
                                   const MAST::FieldFunction<Real>& rho,
                                   const MAST::FieldFunction<Real>& gamma,
                                   const MAST::FieldFunction<Real>& dwdx,
                                   const MAST::FieldFunction<Real>& dwdt):
MAST::FieldFunction<Real>("pressure"),
_order(order),
_V_inf(V),
_M_inf(M),
_rho_inf(rho),
_gamma(gamma),
_dwdx(dwdx),
_dwdt(dwdt) {
    
    _functions.insert(&_V_inf);
    _functions.insert(&_M_inf);
    _functions.insert(&_rho_inf);
    _functions.insert(&_gamma);
    _functions.insert(&_dwdx);
    _functions.insert(&_dwdt);
}




MAST::PistonTheoryPressureXdotDerivative::~PistonTheoryPressureXdotDerivative() { }



void
MAST::PistonTheoryPressureXdotDerivative::operator() (const libMesh::Point& p,
                                                      const Real t,
                                                      Real& m) const {
    
    Real V, M, rho, gamma, dwdx, dwdt;
    _V_inf    (p, t,     V);
    _M_inf    (p, t,     M);
    _rho_inf  (p, t,   rho);
    _gamma    (p, t, gamma);
    _dwdx     (p, t,  dwdx);
    _dwdt     (p, t,  dwdt);
    
    
    Real
    Mf      = (M*M-2)/(M*M-1),
    f1      = rho/sqrt(M*M-1.),

    v0      = dwdx + (M*M-2)/(M*M-1)*dwdt/V,
    v2      = pow(V*dwdx + (M*M-2)/(M*M-1)*dwdt, 2),
    
    v0dwt   = Mf/V,
    v1dwt   = Mf*V,
    v2dwt   = 2.*(V*dwdx + Mf*dwdt)*Mf,
    v3dwt   = v0dwt*v2+v0*v2dwt;
    
    //
    //  linear part:  f1 * dv1/dwt,
    //  where    f1 = rho/sqrt(M*M-1.) and
    //           v1 = V*V*dwdx + Mf*V*dwdt,
    //
    m = f1*v1dwt; // linear
    
    switch (_order) {
            
        case 3:
            //
            //  cubic part:  f1 * (gamma+1)/12.*M*M*dv3/dwt
            //
            m += f1*M*M*v3dwt * (gamma+1)/12.;
        case 2:
            //
            //  quadratic part:  f1 * (gamma+1)/4.*M*dv2/dwt
            //
            m += f1*M*v2dwt * (gamma+1)/4.;
    }
}




void
MAST::PistonTheoryPressureXdotDerivative::derivative (                     const MAST::FunctionBase& f,
                                                      const libMesh::Point& p,
                                                      const Real t,
                                                      Real& m) const {
    
    Real V, M, rho, gamma, dwdx, dwdt, dV, dM, drho, dgamma;
    _V_inf   (p, t,     V);  _V_inf.derivative   (f, p, t,     dV);
    _M_inf   (p, t,     M);  _M_inf.derivative   (f, p, t,     dM);
    _rho_inf (p, t,   rho);  _rho_inf.derivative (f, p, t,   drho);
    _gamma   (p, t, gamma);  _gamma.derivative   (f, p, t, dgamma);
    _dwdx    (p, t,  dwdx);
    _dwdt    (p, t,  dwdt);
    
    Real
    Mf        = (M*M-2)/(M*M-1),
    dMfdM     = 2*M*dM/(M*M-1) - (M*M-2)/pow(M*M-1,2)*2*M*dM,
    
    v0        = dwdx + Mf*dwdt/V,
    v2        = pow(V*dwdx + Mf*dwdt, 2),
    
    dv0dV     = -Mf*dwdt/(V*V)*dV,
    dv2dV     = 2.*(V*dwdx + Mf*dwdt)*dwdx*dV,
    
    dv0dM     = dMfdM*dwdt/V,
    dv2dM     = 2.*(V*dwdx + Mf*dwdt)*dMfdM*dwdt,
    
    v0dwt     = Mf/V,
    v1dwt     = Mf*V,
    v2dwt     = 2.*(V*dwdx + Mf*dwdt)*Mf,
    v3dwt     = v0dwt*v2+v0*v2dwt,
    
    v0dwtdV   = -Mf/(V*V)*dV,
    v1dwtdV   = Mf*dV,
    v2dwtdV   = 2.*dV*dwdx *Mf,
    v3dwtdV   = v0dwtdV*v2+v0dwt*dv2dV + dv0dV*v2dwt+v0*v2dwtdV,
    
    v0dwtdM   = dMfdM/V,
    v1dwtdM   = dMfdM*V,
    v2dwtdM   = 2.*(V*dMfdM*dwdx + 2*Mf*dMfdM*dwdt),
    v3dwtdM   = v0dwtdM*v2+v0dwt*dv2dM + dv0dM*v2dwt+v0*v2dwtdM,

    f1        = rho/sqrt(M*M-1.),
    df1drho   = drho/sqrt(M*M-1.),
    df1dM     = -.5*rho/pow(M*M-1., 1.5)*2*M*dM;
    
    
    //
    //  linear part:  d/dp (f1 * dv1/dwt),
    //  where    f1 = rho/sqrt(M*M-1.) and
    //           v1 = V*V*dwdx + Mf*V*dwdt,
    //
    m = f1*(v1dwtdV + v1dwtdM) + (df1drho + df1dM) * v1dwt; // linear
    
    switch (_order) {
            
        case 3:
            //
            //  cubic part:  d/dp (f1 * (gamma+1)/12.*M*M* dv3/dwt)
            //
            m +=
            (df1drho + df1dM) *M*M*v3dwt * (gamma+1)/12.+
            f1*M*(2*dM*v3dwt + M*(v3dwtdV + v3dwtdM)) * (gamma+1)/12. +
            f1 * M*M*v3dwt/12.*dgamma;
        case 2:
            //
            //  quadratic part:  d/dp (f1 * (gamma+1)/4.*M* dv2/dwt)
            //
            m +=
            (df1drho + df1dM) * M*v2dwt * (gamma+1)/4. +
            f1*(dM*v2dwt + M*(v2dwtdV + v2dwtdM)) * (gamma+1)/4. +
            f1*M*v2dwt/4.*dgamma;
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





std::unique_ptr<MAST::FieldFunction<Real> >
MAST::PistonTheoryBoundaryCondition::
get_pressure_function(const MAST::FieldFunction<Real>& dwdx,
                      const MAST::FieldFunction<Real>& dwdt) const {
    
    MAST::PistonTheoryPressure
    *rval = new MAST::PistonTheoryPressure
    (_order,
     this->get<MAST::FieldFunction<Real> >("V"),
     this->get<MAST::FieldFunction<Real> >("mach"),
     this->get<MAST::FieldFunction<Real> >("rho"),
     this->get<MAST::FieldFunction<Real> >("gamma"),
     dwdx,
     dwdt);
    
    return std::unique_ptr<MAST::FieldFunction<Real> >(rval);
}




std::unique_ptr<MAST::FieldFunction<Real> >
MAST::PistonTheoryBoundaryCondition::
get_dpdx_function(const MAST::FieldFunction<Real>& dwdx,
                  const MAST::FieldFunction<Real>& dwdt) const {
    
    MAST::PistonTheoryPressureXDerivative
    *rval = new MAST::PistonTheoryPressureXDerivative
    (_order,
     this->get<MAST::FieldFunction<Real> >("V"),
     this->get<MAST::FieldFunction<Real> >("mach"),
     this->get<MAST::FieldFunction<Real> >("rho"),
     this->get<MAST::FieldFunction<Real> >("gamma"),
     dwdx,
     dwdt);
    
    return std::unique_ptr<MAST::FieldFunction<Real> >(rval);
}




std::unique_ptr<MAST::FieldFunction<Real> >
MAST::PistonTheoryBoundaryCondition::
get_dpdxdot_function(const MAST::FieldFunction<Real>& dwdx,
                     const MAST::FieldFunction<Real>& dwdt) const {
    
    MAST::PistonTheoryPressureXdotDerivative
    *rval = new MAST::PistonTheoryPressureXdotDerivative
    (_order,
     this->get<MAST::FieldFunction<Real> >("V"),
     this->get<MAST::FieldFunction<Real> >("mach"),
     this->get<MAST::FieldFunction<Real> >("rho"),
     this->get<MAST::FieldFunction<Real> >("gamma"),
     dwdx,
     dwdt);
    
    return std::unique_ptr<MAST::FieldFunction<Real> >(rval);
}


