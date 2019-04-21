/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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

#ifndef __mast_complex_assembly_elem_operations_h__
#define __mast_complex_assembly_elem_operations_h__

// MAST includes
#include "base/assembly_elem_operation.h"
#include "base/mast_data_types.h"

namespace MAST {
    
    class ComplexAssemblyElemOperations:
    public MAST::AssemblyElemOperations {
        
    public:
        ComplexAssemblyElemOperations();
        
        virtual ~ComplexAssemblyElemOperations();

        /*!
         *   sets the element complex solution
         */
        virtual void set_elem_complex_solution(const ComplexVectorX& sol);

        /*!
         *   sets the element complex solution
         */
        virtual void set_elem_complex_solution_sensitivity(const ComplexVectorX& sol);


        /*!
         *   performs the element calculations over \p elem, and returns
         *   the element vector and matrix quantities in \p mat and
         *   \p vec, respectively. \p if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void elem_calculations(bool if_jac,
                                       ComplexVectorX& vec,
                                       ComplexMatrixX& mat) = 0;
        
        /*!
         *   performs the element sensitivity calculations over \p elem,
         *   and returns the element residual sensitivity in \p vec .
         */
        virtual void elem_sensitivity_calculations(const MAST::FunctionBase& f,
                                                   ComplexVectorX& vec) = 0;

    protected:
        
    };

}
#endif // __mast_complex_assembly_elem_operations_h__
