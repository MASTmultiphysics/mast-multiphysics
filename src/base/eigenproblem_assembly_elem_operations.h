/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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

#ifndef __mast_eigenproblem_assembly_elem_operations_h__
#define __mast_eigenproblem_assembly_elem_operations_h__

// MAST includes
#include "base/assembly_elem_operation.h"
#include "base/mast_data_types.h"

namespace MAST {
    
    class EigenproblemAssemblyElemOperations:
    public MAST::AssemblyElemOperations {
        
    public:
        EigenproblemAssemblyElemOperations();
        
        virtual ~EigenproblemAssemblyElemOperations();
        
        /*!
         *   sets the element solution(s) before calculations
         */
        virtual void
        set_elem_sol(MAST::ElementBase& elem,
                     const RealVectorX& sol) = 0;
        
        /*!
         *   sets the element solution(s) before calculations
         */
        virtual void
        set_elem_sol_sens(MAST::ElementBase& elem,
                          const RealVectorX& sol) = 0;
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        elem_calculations(MAST::ElementBase& elem,
                          RealMatrixX& mat_A,
                          RealMatrixX& mat_B) = 0;
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        elem_sensitivity_calculations(MAST::ElementBase& elem,
                                      bool base_sol,
                                      RealMatrixX& mat_A,
                                      RealMatrixX& mat_B) = 0;

    protected:
        
    };
}



#endif // __mast_eigenproblem_assembly_elem_operations_h__
