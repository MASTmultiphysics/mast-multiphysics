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


#ifndef __mast__heat_conduction_nonlinear_assembly__
#define __mast__heat_conduction_nonlinear_assembly__

// MAST includes
#include "base/nonlinear_implicit_assembly_elem_operations.h"


namespace MAST {
    
    
    class HeatConductionNonlinearAssemblyElemOperations:
    public MAST::NonlinearImplicitAssemblyElemOperations {
        
    public:
        
        /*!
         *   constructor associates this assembly object with the system
         */
        HeatConductionNonlinearAssemblyElemOperations();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~HeatConductionNonlinearAssemblyElemOperations();
        
        /*!
         *   performs the element calculations over \p elem, and returns
         *   the element vector and matrix quantities in \p mat and
         *   \p vec, respectively. \p if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void
        elem_calculations(bool if_jac,
                          RealVectorX& vec,
                          RealMatrixX& mat);
        
        /*!
         *   performs the element sensitivity calculations over \p elem,
         *   and returns the element residual sensitivity in \p vec .
         */
        virtual void
        elem_sensitivity_calculations(const MAST::FunctionBase& f,
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
                                               const MAST::LevelSetIntersection& intersect,
                                               const MAST::FieldFunction<RealVectorX>& vel,
                                               RealVectorX& vec);

        /*!
         *   calculates \f$ d ([J] \{\Delta X\})/ dX  \f$ over \p elem,
         *   and returns the matrix in \p vec .
         */
        virtual void
        elem_second_derivative_dot_solution_assembly(RealMatrixX& mat);

        virtual void
        elem_linearized_jacobian_solution_product(RealVectorX& vec) {
            
            libmesh_assert(false); // not implemented yet.
        }

        /*!
         *   initializes the object for the geometric element \p elem. This
         *   expects the object to be in a cleared state, so the user should
         *   call \p clear_elem() between successive initializations.
         */
        virtual void
        init(const libMesh::Elem& elem);
        
    protected:
        
    };
}


#endif // __mast__heat_conduction_nonlinear_assembly__
