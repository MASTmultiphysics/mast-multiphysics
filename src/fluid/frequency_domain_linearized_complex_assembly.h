/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef __mast__frequency_domain_linearized_complex_assembly_elem_operations_h__
#define __mast__frequency_domain_linearized_complex_assembly_elem_operations_h__

// MAST includes
#include "base/complex_assembly_elem_operations.h"



namespace MAST {
    
    
    // Forward declerations
    class FunctionBase;
    class FrequencyFunction;
    
    
    class FrequencyDomainLinearizedComplexAssemblyElemOperations:
    public MAST::ComplexAssemblyElemOperations {
        
    public:
        
        /*!
         *   default constructor
         */
        FrequencyDomainLinearizedComplexAssemblyElemOperations();
        
        
        /*!
         *   destructor
         */
        virtual ~FrequencyDomainLinearizedComplexAssemblyElemOperations();
        

        /*!
         *    sets the frequency function for analysis
         */
        void set_frequency_function(MAST::FrequencyFunction& f);
        
        
        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void clear_frequency_function();

        /*!
         *   performs the element calculations over \p elem, and returns
         *   the element vector and matrix quantities in \p mat and
         *   \p vec, respectively. \p if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void elem_calculations(bool if_jac,
                                       ComplexVectorX& vec,
                                       ComplexMatrixX& mat);
        
        /*!
         *   performs the element sensitivity calculations over \p elem,
         *   and returns the element residual sensitivity in \p vec .
         */
        virtual void elem_sensitivity_calculations(const MAST::FunctionBase& f,
                                                   ComplexVectorX& vec);
        
        /*!
         *   virtual function, but nothing to be done for fluids.
         */
        virtual void
        set_elem_data(unsigned int dim,
                      const libMesh::Elem& ref_elem,
                      MAST::GeomElem& elem) const { }

        /*!
         *   initializes the object for the geometric element \p elem. This
         *   expects the object to be in a cleared state, so the user should
         *   call \p clear_elem() between successive initializations.
         */
        virtual void
        init(const MAST::GeomElem& elem);

    protected:
        
        /*!
         *   frequency function used to define the oscillatory frequency
         */
        MAST::FrequencyFunction*  _frequency;
        
    };
}


#endif // __mast__frequency_domain_linearized_complex_assembly_elem_operations_h__

