
// MAST includes
#include "heat_conduction/heat_conduction_system_initialization.h"


MAST::HeatConductionSystemInitialization::
HeatConductionSystemInitialization(libMesh::System& sys,
                                   const std::string& prefix,
                                   const libMesh::FEType& fe_type):
MAST::SystemInitialization(sys, prefix) {

    _vars.resize(1);
    
    std::string nm = prefix + "_T";
    
    _vars[0] = sys.add_variable(nm, fe_type);
    
}


MAST::HeatConductionSystemInitialization::~HeatConductionSystemInitialization() {
    
}


