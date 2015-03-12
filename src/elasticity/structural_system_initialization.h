
#ifndef __mast__structural_system_initialization__
#define __mast__structural_system_initialization__

// MAST includes
#include "base/system_initialization.h"


namespace MAST {
    
    class StructuralSystemInitialization:
    public MAST::SystemInitialization  {
        
    public:
        StructuralSystemInitialization(libMesh::System& sys,
                                       const std::string& prefix,
                                       const libMesh::FEType& fe_type);
        
        virtual ~StructuralSystemInitialization();
        
    protected:
        
    };
}

#endif // __mast__structural_system_initialization__
