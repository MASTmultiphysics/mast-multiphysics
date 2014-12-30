
// MAST includes
#include "base/assembly_base.h"
#include "base/system_initialization.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"


MAST::AssemblyBase::AssemblyBase():
_discipline(NULL),
_system(NULL),
_sol_function(NULL) {
    
}




MAST::AssemblyBase::~AssemblyBase() {
    
}



const MAST::PhysicsDisciplineBase&
MAST::AssemblyBase::discipline() const {
    
    libmesh_assert_msg(_discipline,
                       "Error: Discipline not yet attached to Assembly.");
    return *_discipline;
}



MAST::PhysicsDisciplineBase&
MAST::AssemblyBase::discipline() {
    
    libmesh_assert_msg(_discipline,
                       "Error: Discipline not yet attached to Assembly.");
    return *_discipline;
}




const libMesh::System&
MAST::AssemblyBase::system() const {
    
    libmesh_assert_msg(_discipline,
                       "Error: System not yet attached to Assembly.");
    return _system->system();
}


libMesh::System&
MAST::AssemblyBase::system() {
    
    libmesh_assert_msg(_discipline,
                       "Error: System not yet attached to Assembly.");
    return _system->system();
}




std::auto_ptr<libMesh::NumericVector<Real> >
MAST::AssemblyBase::_build_localized_vector(const libMesh::System& sys,
                                            const libMesh::NumericVector<Real>& global) {
    
    libMesh::NumericVector<Real>* local =
    libMesh::NumericVector<Real>::build(sys.comm()).release();
    
    const std::vector<libMesh::dof_id_type>& send_list =
    sys.get_dof_map().get_send_list();
    
    local->init(sys.n_dofs(),
                sys.n_local_dofs(),
                send_list,
                false,
                libMesh::GHOSTED);
    global.localize(*local, send_list);
    
    return std::auto_ptr<libMesh::NumericVector<Real> >(local);
}




void
MAST::AssemblyBase::attach_solution_function(MAST::MeshFieldFunction<RealVectorX>& f){
    
    // make sure that no prior association is specified
    libmesh_assert(!_sol_function);
    
    _sol_function = &f;
}




void
MAST::AssemblyBase::detach_solution_function() {
    _sol_function = NULL;
}


