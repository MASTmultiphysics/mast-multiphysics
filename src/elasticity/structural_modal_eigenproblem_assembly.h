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
        set_elem_solution(const RealVectorX& sol);

        /*!
         *   sets the element solution sensitivity before calculations
         */
        virtual void
        set_elem_solution_sensitivity(const RealVectorX& sol);

        /*!
         *   performs the element calculations over \p elem, and returns
         *   the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        elem_calculations(RealMatrixX& mat_A,
                          RealMatrixX& mat_B);
        
        /*!
         *   performs the element sensitivity calculations over \p elem,
         *   and returns the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        elem_sensitivity_calculations(const MAST::FunctionBase& f,
                                      bool base_sol,
                                      RealMatrixX& mat_A,
                                      RealMatrixX& mat_B);

        /*!
         *   performs the element topology sensitivity calculations over \p elem.
         */
        virtual void
        elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                               bool base_sol,
                                               RealMatrixX& mat_A,
                                               RealMatrixX& mat_B);


        /*!
         *   performs the element topology sensitivity calculations over \p elem.
         */
        virtual void
        elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                               bool base_sol,
                                               const MAST::FieldFunction<RealVectorX>& vel,
                                               RealMatrixX& mat_A,
                                               RealMatrixX& mat_B);

        /*!
         *   sets the structural element y-vector if 1D element is used.
         */
        virtual void
        set_elem_data(unsigned int dim,
                      const libMesh::Elem& ref_elem,
                      MAST::GeomElem& elem) const;

        /*!
         *   initializes the object for the geometric element \p elem. This
         *   expects the object to be in a cleared state, so the user should
         *   call \p clear_elem() between successive initializations.
         */
        virtual void
        init(const MAST::GeomElem& elem);

    protected:
        
        
        MAST::StructuralAssembly* _incompatible_sol_assembly;
    };
    
}

#endif // __mast__structural_modal_eigenproblem_assembly_elem_operations_h__
