
// MAST includes
#include "coordinates/basis_matrix_coordinate.h"


MAST::BasisMatrixCoordinate::BasisMatrixCoordinate(const std::string& nm,
                                                   MAST::FieldFunction<RealMatrixX>* basis):
MAST::CoordinateBase(nm),
_basis(basis) {
    
    _functions.insert(basis->master());
}




MAST::BasisMatrixCoordinate::~BasisMatrixCoordinate() {
    
    delete _basis;
}




MAST::BasisMatrixCoordinate::BasisMatrixCoordinate(const MAST::BasisMatrixCoordinate& c):
MAST::CoordinateBase(c),
_basis(c._basis->clone().release()) {
    
    _functions.insert(_basis->master());
}




std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::BasisMatrixCoordinate::clone() const {
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
    (new MAST::BasisMatrixCoordinate(*this));
}




void
MAST::BasisMatrixCoordinate::operator() (const libMesh::Point& p,
                                         const Real t,
                                         RealMatrixX& v) const {
    
    (*_basis)(p, t, v);
}




void
MAST::BasisMatrixCoordinate::derivative (const MAST::DerivativeType d,
                                         const MAST::FunctionBase& f,
                                         const libMesh::Point& p,
                                         const Real t,
                                         RealMatrixX& v) const {
    
    _basis->derivative(d, f, p, t, v);
}



