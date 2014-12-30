

// MAST includes
#include "mesh/local_elem_base.h"
#include "base/field_function_base.h"


MAST::LocalElemBase::~LocalElemBase() {

}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::LocalElemBase::T_matrix_function() const {

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::TransformMatrixFunction(_T_mat);
        
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}




void
MAST::LocalElemBase::global_coordinates_location(const libMesh::Point& local,
                                                 libMesh::Point& global) const {
    if (!_local_elem) // no local elem is created for the case of 3D elem
        global = local;
    else {
        global = 0.;
        
        // now calculate the global coordinates with respect to the origin
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++)
                global(j) += _T_mat(j,k)*local(k);
        
        // shift to the global coordinate
        global += (*_elem.get_node(0));
    }
}



void
MAST::LocalElemBase::global_coordinates_normal(const libMesh::Point& local,
                                               RealVector3& global) const {
    if (!_local_elem)
        for (unsigned int i=0; i<3; i++)
            global(i) = local(i);
    else {
        global.setZero();
        
        // now calculate the global coordinates with respect to the origin
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++)
                global(j) += _T_mat(j,k)*local(k);
    }
}


MAST::TransformMatrixFunction::
TransformMatrixFunction(const RealMatrixX& Tmat):
MAST::FieldFunction<RealMatrixX>("T_function"),
_Tmat(Tmat) {
    
}



MAST::TransformMatrixFunction::
TransformMatrixFunction(const MAST::TransformMatrixFunction& f):
MAST::FieldFunction<RealMatrixX>(f),
_Tmat(f._Tmat) {
    
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::TransformMatrixFunction::clone() const {
    MAST::FieldFunction<RealMatrixX> *rval =
    new MAST::TransformMatrixFunction(*this);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}



void
MAST::TransformMatrixFunction::operator() (const libMesh::Point& p,
                                           const Real t,
                                           RealMatrixX& v) const {
    v = _Tmat;
}



void
MAST::TransformMatrixFunction::derivative (const MAST::DerivativeType d,
                                           const MAST::FunctionBase& f,
                                           const libMesh::Point& p,
                                           const Real t,
                                           RealMatrixX& v) const {
    v.setZero(3, 3);
}


