

// MAST includes
#include "elasticity/structural_system_initialization.h"


MAST::StructuralSystemInitialization::
StructuralSystemInitialization(libMesh::System& sys,
                               const std::string& prefix,
                               const libMesh::FEType& fe_type):
MAST::SystemInitialization(sys, prefix) {
    
    _vars.resize(6);
    
    std::string nm = prefix + "_ux";
    _vars[0] = sys.add_variable(nm, fe_type);

    nm = prefix + "_uy";
    _vars[1] = sys.add_variable(nm, fe_type);

    nm = prefix + "_uz";
    _vars[2] = sys.add_variable(nm, fe_type);

    nm = prefix + "_tx";
    _vars[3] = sys.add_variable(nm, fe_type);

    nm = prefix + "_ty";
    _vars[4] = sys.add_variable(nm, fe_type);

    nm = prefix + "_tz";
    _vars[5] = sys.add_variable(nm, fe_type);

}


MAST::StructuralSystemInitialization::~StructuralSystemInitialization() {
    
}
