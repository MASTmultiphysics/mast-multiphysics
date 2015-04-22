
#ifndef __mast__structural_modal_eigenproblem_assembly__
#define __mast__structural_modal_eigenproblem_assembly__

// MAST includes
#include "base/eigenproblem_assembly.h"


namespace MAST {
    
    
    class StructuralModalEigenproblemAssembly:
    public MAST::EigenproblemAssembly {
    public:
        
        /*!
         *   constructor associates the eigen system with this assembly object
         */
        StructuralModalEigenproblemAssembly();
        
        /*!
         *   destructor resets the association with the eigen system
         *   from this assembly object
         */
        virtual ~StructuralModalEigenproblemAssembly();

        
        /*!
         *   assembles the global A and B matrices for the modal 
         *   eigenvalue problem
         */
        virtual void assemble();
        
    protected:
        
        /*!
         *   @returns a smart-pointer to a newly created element for
         *   calculation of element quantities.
         */
        virtual std::auto_ptr<MAST::ElementBase>
        _build_elem(const libMesh::Elem& elem);

        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        _elem_calculations(MAST::ElementBase& elem,
                           RealMatrixX& mat_A,
                           RealMatrixX& mat_B);
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                       RealMatrixX& mat_A,
                                       RealMatrixX& mat_B);
        
    };
    
}

#endif // __mast__structural_modal_eigenproblem_assembly__
