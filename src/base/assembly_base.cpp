
// MAST includes
#include "base/assembly_base.h"
#include "base/system_initialization.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"


MAST::AssemblyBase::AssemblyBase(MAST::PhysicsDisciplineBase& discipline,
                                 MAST::SystemInitialization& system):
libMesh::System::QOI(),
libMesh::System::QOIDerivative(),
libMesh::System::QOIParameterSensitivity(),
_discipline(discipline),
_system(system) {
    
    libMesh::System& sys = system.system();
    sys.attach_QOI_object(*this);
    sys.attach_QOI_derivative_object(*this);
    sys.attach_QOI_parameter_sensitivity_object(*this);
}




MAST::AssemblyBase::~AssemblyBase() {
    
    libMesh::System& sys = _system.system();
    
    sys.reset_QOI();
    sys.reset_QOI_derivative();
    sys.reset_QOI_parameter_sensitivity();
}



const libMesh::System&
MAST::AssemblyBase::system() const {
    
    return _system.system();
}


libMesh::System&
MAST::AssemblyBase::system() {
    
    return _system.system();
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





