
#ifndef __mast__heat_conduction_transient_assembly__
#define __mast__heat_conduction_transient_assembly__

// MAST includes
#include "base/transient_assembly.h"



namespace MAST {
    
    
    class HeatConductionTransientAssembly:
    public MAST::TransientAssembly {
    public:
        
        /*!
         *   constructor associates this assembly object with the system
         */
        HeatConductionTransientAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~HeatConductionTransientAssembly();
        
        //**************************************************************
        //these methods are provided for use by the solvers
        //**************************************************************
        
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void _elem_calculations(MAST::ElementBase& elem,
                                        bool if_jac,
                                        RealVectorX& f_m,
                                        RealVectorX& f_x,
                                        RealMatrixX& f_m_jac_x_dot,
                                        RealMatrixX& f_m_jac,
                                        RealMatrixX& f_x_jac);
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                                    RealVectorX& vec);
        
        
    protected:
        
        
        /*!
         *   @returns a smart-pointer to a newly created element for
         *   calculation of element quantities.
         */
        virtual std::auto_ptr<MAST::ElementBase>
        _build_elem(const libMesh::Elem& elem);

    };
    
    
}

#endif // __mast__heat_conduction_transient_assembly__
