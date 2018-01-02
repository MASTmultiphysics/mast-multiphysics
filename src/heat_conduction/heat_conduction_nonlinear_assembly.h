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
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void
        elem_calculations(MAST::ElementBase& elem,
                          bool if_jac,
                          RealVectorX& vec,
                          RealMatrixX& mat);
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void
        elem_sensitivity_calculations(MAST::ElementBase& elem,
                                      RealVectorX& vec);
        
        /*!
         *   calculates \f$ d ([J] \{\Delta X\})/ dX  \f$ over \par elem,
         *   and returns the matrix in \par vec .
         */
        virtual void
        elem_second_derivative_dot_solution_assembly(MAST::ElementBase& elem,
                                                      RealMatrixX& mat);

        virtual void
        elem_linearized_jacobian_solution_product(MAST::ElementBase& elem,
                                                  RealVectorX& vec) {
            libmesh_assert(false); // not implemented yet.
        }

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
        
    };
}


#endif // __mast__heat_conduction_nonlinear_assembly__
