
#ifndef __mast__driver_base__
#define __mast__driver_base__

namespace MAST {
    // Forward declerations
    class PhysicsDisciplineBase;
    class SystemInitialization;
    class NonlinearImplicitAssembly;
    class OutputAssemblyBase;
    

    namespace Driver {
        
        bool nonlinear_solution(MAST::PhysicsDisciplineBase&     discipline,
                                MAST::SystemInitialization&      system,
                                MAST::NonlinearImplicitAssembly& assembly);

        bool output_evaluation(MAST::PhysicsDisciplineBase&     discipline,
                               MAST::SystemInitialization&      system,
                               MAST::OutputAssemblyBase&        assembly);

    
        bool adjoint_solution(MAST::PhysicsDisciplineBase&     discipline,
                              MAST::SystemInitialization&      system,
                              MAST::NonlinearImplicitAssembly& nonlin_assembly,
                              MAST::OutputAssemblyBase&        output_assembly);
    }
}


#endif // __mast__driver_base__
