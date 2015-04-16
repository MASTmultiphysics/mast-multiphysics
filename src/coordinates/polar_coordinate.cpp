
// MAST includes
#include "coordinates/polar_coordinate.h"


MAST::PolarCoordinate::PolarCoordinate(const std::string& nm,
                                       MAST::FieldFunction<Real>* theta):
MAST::CoordinateBase(nm),
_theta(theta) {
    
    _functions.insert(theta->master());
}


MAST::PolarCoordinate::~PolarCoordinate() {
    delete _theta;
}



MAST::PolarCoordinate::PolarCoordinate(const MAST::PolarCoordinate& c):
MAST::CoordinateBase(c),
_theta(c._theta->clone().release()) {
    _functions.insert(_theta->master());
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::PolarCoordinate::clone() const {
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
    (new MAST::PolarCoordinate(*this));
}



void
MAST::PolarCoordinate::operator() (const libMesh::Point& p,
                                   const Real t,
                                   RealMatrixX& v) const {
    v.setZero(3,3);
    Real theta;
    (*_theta)  (p,t,theta);
    v(0,0) = cos(theta);
    v(1,0) = sin(theta);
    
    v(0,1) = -sin(theta);
    v(1,1) =  cos(theta);
    
    v(2,2) = 1.;
}



void
MAST::PolarCoordinate::derivative (const MAST::DerivativeType d,
                                   const MAST::FunctionBase& f,
                                   const libMesh::Point& p,
                                   const Real t,
                                   RealMatrixX& v) const {
    v.setZero(3,3);
    Real theta, dtheta;
    (*_theta)  (p,t,theta);
    _theta->derivative(d,f,p,t,dtheta);
    
    v(0,0) = -sin(theta)*dtheta;
    v(1,0) =  cos(theta)*dtheta;
    
    v(0,1) = -cos(theta)*dtheta;
    v(1,1) = -sin(theta)*dtheta;
}
