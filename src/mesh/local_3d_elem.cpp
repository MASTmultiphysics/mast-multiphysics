

// MAST includes
#include "mesh/local_3d_elem.h"


MAST::Local3DElem::Local3DElem(const libMesh::Elem& elem):
MAST::LocalElemBase(elem) {
    
}



MAST::Local3DElem::~Local3DElem() {
    
}



void
MAST::Local3DElem::
domain_surface_normal_in_global_coordinates(const libMesh::Point& p,
                                            RealVector3& n_global) const {
    for (unsigned int i=0; i<3; i++)
        n_global(i) = p(i);
}


