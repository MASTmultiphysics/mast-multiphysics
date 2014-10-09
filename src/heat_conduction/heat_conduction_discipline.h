
#ifndef __mast__heat_conduction_discipline__
#define __mast__heat_conduction_discipline__

// MAST includes
#include "base/physics_discipline_base.h"



namespace MAST {
    
    
    class HeatConductionDiscipline:
    public MAST::PhysicsDisciplineBase {
        

    public:
        
        // Constructor
        HeatConductionDiscipline(libMesh::EquationSystems& eq_sys);
        
        
        /*!
         *   virtual destructor
         */
        virtual ~HeatConductionDiscipline();

        
    protected:
        
        
    };
}


#endif // __mast__heat_conduction_discipline__
