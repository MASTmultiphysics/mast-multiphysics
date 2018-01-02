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

#ifndef __mast__structural_modal_eigenproblem_assembly_elem_operations_h__
#define __mast__structural_modal_eigenproblem_assembly_elem_operations_h__

// MAST includes
#include "base/eigenproblem_assembly_elem_operations.h"


namespace MAST {
    
    // Forward declerations
    class StructuralAssembly;
    
    class StructuralModalEigenproblemAssemblyElemOperations:
    public MAST::EigenproblemAssemblyElemOperations {
    public:
        
        /*!
         *   constructor associates the eigen system with this assembly object
         */
        StructuralModalEigenproblemAssemblyElemOperations();
        
        /*!
         *   destructor resets the association with the eigen system
         *   from this assembly object
         */
        virtual ~StructuralModalEigenproblemAssemblyElemOperations();

                
        /*!
         *   attached the incompatible solution object
         */
        void attach_incompatible_solution_object(MAST::StructuralAssembly& str_assembly);
        
        
        /*!
         *   clear the incompatible solution object
         */
        void clear_incompatible_solution_object();

        /*!
         *   sets the element solution(s) before calculations
         */
        virtual void
        set_elem_sol(MAST::ElementBase& elem,
                     const RealVectorX& sol);

        /*!
         *   sets the element solution sensitivity before calculations
         */
        virtual void
        set_elem_sol_sens(MAST::ElementBase& elem,
                          const RealVectorX& sol);

        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        elem_calculations(MAST::ElementBase& elem,
                          RealMatrixX& mat_A,
                          RealMatrixX& mat_B);
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        elem_sensitivity_calculations(MAST::ElementBase& elem,
                                      bool base_sol,
                                      RealMatrixX& mat_A,
                                      RealMatrixX& mat_B);

        
        /*!
         *   @returns a MAST::FEBase object for calculation of finite element
         *   quantities. This creates LocalElemFE for 1D and 2D elements.
         */
        virtual std::unique_ptr<MAST::FEBase>
        build_fe(const libMesh::Elem& e);

        /*!
         *   @returns a smart-pointer to a newly created element for
         *   calculation of element quantities.
         */
        virtual std::unique_ptr<MAST::ElementBase>
        build_elem(const libMesh::Elem& elem);

    protected:
        
        
        MAST::StructuralAssembly* _incompatible_sol_assembly;
    };
    
}

#endif // __mast__structural_modal_eigenproblem_assembly_elem_operations_h__
