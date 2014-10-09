
// MAST includes
#include "base/output_assembly_base.h"
#include "base/system_initialization.h"




MAST::OutputAssemblyBase::OutputAssemblyBase():
MAST::AssemblyBase(),
libMesh::System::QOI(),
libMesh::System::QOIDerivative(),
libMesh::System::QOIParameterSensitivity() {
    
}



MAST::OutputAssemblyBase::~OutputAssemblyBase() {
    
}



void
MAST::OutputAssemblyBase::
attach_discipline_and_system(MAST::PhysicsDisciplineBase &discipline,
                             MAST::SystemInitialization &system) {
    
    libmesh_assert_msg(!_discipline && !_system,
                       "Error: Assembly should be cleared before attaching System.");
    
    _discipline = &discipline;
    _system     = &system;
    
    // now attach this to the system
    libMesh::System& sys = system.system();
    sys.attach_QOI_object(*this);
    sys.attach_QOI_derivative_object(*this);
    sys.attach_QOI_parameter_sensitivity_object(*this);
}




void
MAST::OutputAssemblyBase::
clear_discipline_and_system() {
    
    if (_system && _discipline) {
        libMesh::System& sys = _system->system();
        
        sys.reset_QOI();
        sys.reset_QOI_derivative();
        sys.reset_QOI_parameter_sensitivity();
    }
    
    _discipline = NULL;
    _system     = NULL;
}





