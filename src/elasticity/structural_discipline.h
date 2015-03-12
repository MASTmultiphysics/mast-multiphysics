
#ifndef __mast__structural_discipline__
#define __mast__structural_discipline__

// MAST includes
#include "base/physics_discipline_base.h"



namespace MAST {
    
    
    class StructuralDiscipline:
    public MAST::PhysicsDisciplineBase {
        
        
    public:
        
        // Constructor
        StructuralDiscipline(libMesh::EquationSystems& eq_sys);
        
        
        /*!
         *   virtual destructor
         */
        virtual ~StructuralDiscipline();
        
        
    protected:
        
        
    };
}

#endif // __mast__structural_discipline__
