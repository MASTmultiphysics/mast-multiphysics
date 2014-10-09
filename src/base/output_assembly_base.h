
#ifndef __mast__output_assembly_base__
#define __mast__output_assembly_base__

// MAST includes
#include "base/assembly_base.h"


namespace MAST {
    
    
    class OutputAssemblyBase:
    public MAST::AssemblyBase,
    public libMesh::System::QOI,
    public libMesh::System::QOIDerivative,
    public libMesh::System::QOIParameterSensitivity {
    
    public:
        OutputAssemblyBase();
        
        virtual ~OutputAssemblyBase();
        
        /*!
         *   attaches a system to this discipline, and vice-a-versa
         */
        virtual void
        attach_discipline_and_system(MAST::PhysicsDisciplineBase& discipline,
                                     MAST::SystemInitialization& system) ;
        
        
        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void
        clear_discipline_and_system( );
        
        
        
    protected:
        
    };
}


#endif // __mast__output_assembly_base__

