

#ifndef __mast__bending_structural_element__
#define __mast__bending_structural_element__

// MAST includes
#include "elasticity/structural_element_base.h"


namespace MAST {
    
    // Forward declerations
    class BendingOperator;
    class BoundaryConditionBase;
    class FEMOperatorMatrix;
    
    
    
    class BendingStructuralElem:
    public MAST::StructuralElementBase {
        
    public:
        BendingStructuralElem(MAST::SystemInitialization& sys,
                              const libMesh::Elem& elem,
                              const MAST::ElementPropertyCardBase& p);
        
        /*!
         *    row dimension of the direct strain matrix, also used for the
         *    bending operator row dimension
         */
        virtual unsigned int n_direct_strain_components() = 0;
        
        /*!
         *    row dimension of the von Karman strain matrix
         */
        virtual unsigned int n_von_karman_strain_components() = 0;
        
    protected:
        
        
        /*!
         *   initialize membrane strain operator matrix
         */
        virtual void initialize_direct_strain_operator(const unsigned int qp,
                                                       FEMOperatorMatrix& Bmat) = 0;
        
        
        /*!
         *    bending operator used for this elmeent
         */
        std::auto_ptr<MAST::BendingOperator> _bending_operator;
        
    };
}


#endif // __mast__bending_structural_element__
