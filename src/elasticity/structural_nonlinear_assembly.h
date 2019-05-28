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

#ifndef __mast__structural_nonlinear_assembly_elem_operations__
#define __mast__structural_nonlinear_assembly_elem_operations__

// MAST includes
#include "base/nonlinear_implicit_assembly_elem_operations.h"


namespace MAST {
    

    // Forward declerations
    class StructuralAssembly;
    
    
    class StructuralNonlinearAssemblyElemOperations:
    public MAST::NonlinearImplicitAssemblyElemOperations {
        
    public:
        
        /*!
         *   constructor associates this assembly object with the system
         */
        StructuralNonlinearAssemblyElemOperations();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~StructuralNonlinearAssemblyElemOperations();
        
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
        virtual void set_elem_solution(const RealVectorX& sol);
        
        /*!
         *   performs the element calculations over \p elem, and returns
         *   the element vector and matrix quantities in \p mat and
         *   \p vec, respectively. \p if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void elem_calculations(bool if_jac,
                                       RealVectorX& vec,
                                       RealMatrixX& mat);
        
        
        /*!
         *   performs the element calculations over \p elem, and returns
         *   the element vector quantity in \p vec. The vector quantity only
         *   include the \f$ [J] \{dX\} \f$ components, so the inherited classes
         *   must ensure that no component of constant forces (traction/body
         *   forces/etc.) are added to this vector.
         */
        virtual void
        elem_linearized_jacobian_solution_product(RealVectorX& vec);
        
        
        /*!
         *   performs the element sensitivity calculations over \p elem,
         *   and returns the element residual sensitivity in \p vec .
         */
        virtual void elem_sensitivity_calculations(const MAST::FunctionBase& f,
                                                   RealVectorX& vec);
        
        /*!
         *   performs the element shape sensitivity calculations over \p elem,
         *   and returns the element residual sensitivity in \p vec .
         */
        virtual void
        elem_shape_sensitivity_calculations(const MAST::FunctionBase& f,
                                            RealVectorX& vec) {
            libmesh_assert(false); // to be implemented
        }
        
        /*!
         *   performs the element topology sensitivity calculations over \p elem,
         *   and returns the element residual sensitivity in \p vec .
         */
        virtual void
        elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                               const MAST::FieldFunction<RealVectorX>& vel,
                                               RealVectorX& vec);

        /*!
         *   calculates \f$ d ([J] \{\Delta X\})/ dX  \f$ over \p elem,
         *   and returns the matrix in \p vec .
         */
        virtual void
        elem_second_derivative_dot_solution_assembly(RealMatrixX& mat);

        
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


#endif // __mast__structural_nonlinear_assembly_elem_operations__
