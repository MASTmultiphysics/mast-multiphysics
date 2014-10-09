
#ifndef __mast__heat_conduction_system_initialization__
#define __mast__heat_conduction_system_initialization__

// MAST includes
#include "base/system_initialization.h"


namespace MAST {
    
    class HeatConductionSystemInitialization:
    public MAST::SystemInitialization  {
        
    public:
        HeatConductionSystemInitialization(libMesh::System& sys,
                                           const std::string& prefix,
                                           const libMesh::FEType& fe_type);
        
        virtual ~HeatConductionSystemInitialization();
        
    protected:

    };
}

#endif // __mast__heat_conduction_system_initialization__
