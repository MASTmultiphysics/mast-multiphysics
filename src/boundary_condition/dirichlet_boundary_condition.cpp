
// MAST includes
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "base/field_function_base.h"


// libMesh includes
#include "libmesh/zero_function.h"

void
MAST::DirichletBoundaryCondition::init(const libMesh::boundary_id_type bid,
                                       const std::vector<unsigned int>& constrained_vars,
                                       MAST::FieldFunction<Real>* f) {
    
    // should not have been initialized if this is called
    libmesh_assert(_dirichlet_boundary.get() == NULL);
    
    std::set<libMesh::boundary_id_type> bid_set;
    bid_set.insert(bid);
    
    std::auto_ptr<libMesh::FunctionBase<Real> > function;

    // if the function was not give, then assume it to be zero function
    if (f)
        function.reset(f->libMesh_compatible_function().release());
    else
        function.reset(new libMesh::ZeroFunction<Real>);
    
    _dirichlet_boundary.reset(new libMesh::DirichletBoundary(bid_set,
                                                             constrained_vars,
                                                             function.get()));
}


