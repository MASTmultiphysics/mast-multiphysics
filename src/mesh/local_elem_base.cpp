

// MAST includes
#include "mesh/local_elem_base.h"
#include "base/field_function_base.h"


MAST::LocalElemBase::~LocalElemBase() {
    delete _T_mat_function;
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::LocalElemBase::T_matrix_function() const {
    
    return _T_mat_function->clone();
}




void
MAST::LocalElemBase::global_coordinates_location(const libMesh::Point& local,
                                                 libMesh::Point& global) const {
    global = 0.;
    
    // now calculate the global coordinates with respect to the origin
    for (unsigned int j=0; j<3; j++)
        for (unsigned int k=0; k<3; k++)
            global(j) += _T_mat(j,k)*local(k);
    
    // shift to the global coordinate
    global += (*_elem.get_node(0));
}



void
MAST::LocalElemBase::global_coordinates_normal(const libMesh::Point& local,
                                               RealVector3& global) const {
    global.setZero();
    
    // now calculate the global coordinates with respect to the origin
    for (unsigned int j=0; j<3; j++)
        for (unsigned int k=0; k<3; k++)
            global(j) += _T_mat(j,k)*local(k);
}
