
#ifndef __mast__driver_base__
#define __mast__driver_base__

namespace MAST {
    // Forward declerations
    class NonlinearImplicitAssembly;
    
    namespace Driver {
        
        bool nonlinear_solution(MAST::NonlinearImplicitAssembly& assembly);
    }
}


#endif // __mast__driver_base__
