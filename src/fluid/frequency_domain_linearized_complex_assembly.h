/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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

#ifndef __mast__frequency_domain_linearized_complex_assembly_h__
#define __mast__frequency_domain_linearized_complex_assembly_h__

// MAST includes
#include "base/complex_assembly_base.h"



namespace MAST {
    
    
    // Forward declerations
    class FunctionBase;
    class FrequencyFunction;
    
    
    class FrequencyDomainLinearizedComplexAssembly:
    public MAST::ComplexAssemblyBase {
        
    public:
        
        /*!
         *   default constructor
         */
        FrequencyDomainLinearizedComplexAssembly();
        
        
        /*!
         *   destructor
         */
        virtual ~FrequencyDomainLinearizedComplexAssembly();
        

        /*!
         *    sets the frequency function for analysis
         */
        void set_frequency_function(MAST::FrequencyFunction& f);
        
        
        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void clear_discipline_and_system( );

        
    protected:
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void _elem_calculations(MAST::ElementBase& elem,
                                        bool if_jac,
                                        ComplexVectorX& vec,
                                        ComplexMatrixX& mat);
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                                    bool if_jac,
                                                    ComplexVectorX& vec,
                                                    ComplexMatrixX& mat);

        /*!
         *   @returns a smart-pointer to a newly created element for
         *   calculation of element quantities.
         */
        virtual std::unique_ptr<MAST::ElementBase>
        _build_elem(const libMesh::Elem& elem);

        
        /*!
         *   frequency function used to define the oscillatory frequency
         */
        MAST::FrequencyFunction*  _frequency;
        
    };
}


#endif // __mast__frequency_domain_linearized_complex_assembly_h__
