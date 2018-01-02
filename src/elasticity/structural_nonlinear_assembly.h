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
        virtual void set_elem_sol(MAST::ElementBase& elem,
                                  const RealVectorX& sol);
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void elem_calculations(MAST::ElementBase& elem,
                                       bool if_jac,
                                       RealVectorX& vec,
                                       RealMatrixX& mat);
        
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector quantity in \par vec. The vector quantity only
         *   include the \f$ [J] \{dX\} f$ components, so the inherited classes
         *   must ensure that no component of constant forces (traction/body
         *   forces/etc.) are added to this vector.
         */
        virtual void
        elem_linearized_jacobian_solution_product(MAST::ElementBase& elem,
                                                  RealVectorX& vec);
        
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void elem_sensitivity_calculations(MAST::ElementBase& elem,
                                                   RealVectorX& vec);
        
        
        /*!
         *   calculates \f$ d ([J] \{\Delta X\})/ dX  \f$ over \par elem,
         *   and returns the matrix in \par vec .
         */
        virtual void
        elem_second_derivative_dot_solution_assembly(MAST::ElementBase& elem,
                                                     RealMatrixX& mat);

        
        
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

//        /*!
//         *   Evaluates the volume and boundary outputs for the specified
//         *   solution. This reimplements the virtual method from the parent
//         *   class to handle the structural compliance evaluation.
//         */
//        virtual void calculate_outputs(const libMesh::NumericVector<Real>& X);
//
//
//        /*!
//         *   This reimplements the virtual method from the parent
//         *   class to handle the structural compliance evaluation.
//         */
//        void calculate_output_sensitivity(libMesh::ParameterVector& params,
//                                          const bool if_total_sensitivity,
//                                          const libMesh::NumericVector<Real>& X);
//
//
//        /*!
//         *  calculates the elastic compliance of the system \p S about the
//         *  solution defined by \p X. If sensitivity of the quantity is
//         *  desired with respect to a parameter, then the parameter can
//         *  be specified. If \p dX is provided, then the total sensitivity
//         *  is evaluated with respect to the parameter, otherwise the partial
//         *  derivative of the output quantity is evaluated.
//         */
//        virtual void
//        calculate_compliance (const libMesh::NumericVector<Real>& X,
//                              libMesh::NonlinearImplicitSystem& S,
//                              MAST::RealOutputFunction& output,
//                              const MAST::FunctionBase* f = nullptr,
//                              const libMesh::NumericVector<Real>* dX = nullptr);


    protected:
        
        
        MAST::StructuralAssembly* _incompatible_sol_assembly;
        
    };
}


#endif // __mast__structural_nonlinear_assembly_elem_operations__
