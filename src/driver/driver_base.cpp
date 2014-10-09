
// MAST includes
#include "driver/driver_base.h"
#include "base/nonlinear_implicit_assembly.h"
#include "base/physics_discipline_base.h"

// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"



bool
MAST::Driver::nonlinear_solution(MAST::NonlinearImplicitAssembly& assembly) {
    
    libMesh::NonlinearImplicitSystem&      nonlin_sys   =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(assembly.system());
        
    nonlin_sys.solve();
    
    return true;
}


