
// MAST includes
#include "base/constant_field_function.h"
#include "base/parameter.h"



MAST::ConstantFieldFunction::
ConstantFieldFunction(const std::string& nm,
                      const MAST::Parameter& p):
MAST::FieldFunction<Real>(nm),
_p(p) {
    
}


MAST::ConstantFieldFunction::~ConstantFieldFunction() {
    
}



MAST::ConstantFieldFunction::ConstantFieldFunction(const MAST::ConstantFieldFunction& f):
MAST::FieldFunction<Real>(f),
_p(f._p) {
    
}



std::auto_ptr<MAST::FieldFunction<Real> >
MAST::ConstantFieldFunction::clone() const {
    MAST::FieldFunction<Real>* rval =
    new MAST::ConstantFieldFunction(*this);
    
    return std::auto_ptr<MAST::FieldFunction<Real> >(rval);
}



void
MAST::ConstantFieldFunction::operator() (const libMesh::Point& p,
                                         const Real t,
                                         Real& v) const {
    v = _p();
}




void
MAST::ConstantFieldFunction::derivative (const MAST::DerivativeType d,
                                         const MAST::FunctionBase& f,
                                         const libMesh::Point& p,
                                         const Real t,
                                         Real& v) const {
    
    v = _p.depends_on(f)?1:0;
}

