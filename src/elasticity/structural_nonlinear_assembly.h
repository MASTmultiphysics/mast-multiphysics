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

#ifndef __mast__structural_nonlinear_assembly__
#define __mast__structural_nonlinear_assembly__

// MAST includes
#include "base/nonlinear_implicit_assembly.h"


namespace MAST {
    

    // Forward declerations
    class RealOutputFunction;
    class FunctionBase;
    
    
    class StructuralNonlinearAssembly:
    public MAST::NonlinearImplicitAssembly {
        
    public:
        
        /*!
         *   constructor associates this assembly object with the system
         */
        StructuralNonlinearAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~StructuralNonlinearAssembly();
        

        /*!
         *   attaches a system to this discipline, and vice-a-versa
         */
        virtual void
        attach_discipline_and_system(MAST::PhysicsDisciplineBase& discipline,
                                     MAST::SystemInitialization& system);
        
        
        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void
        clear_discipline_and_system( );

        /*!
         *    function that assembles the matrices and vectors quantities for
         *    nonlinear solution
         */
        virtual void
        residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                               libMesh::NumericVector<Real>* R,
                               libMesh::SparseMatrix<Real>*  J,
                               libMesh::NonlinearImplicitSystem& S);
        
        /**
         * Assembly function.  This function will be called
         * to assemble the RHS of the sensitivity equations (which is -1 times
         * sensitivity of system residual) prior to a solve and must
         * be provided by the user in a derived class. The method provides dR/dp_i
         * for \par i ^th parameter in the vector \par parameters.
         *
         * If the routine is not able to provide sensitivity for this parameter,
         * then it should return false, and the system will attempt to use
         * finite differencing.
         */
        virtual bool
        sensitivity_assemble (const libMesh::ParameterVector& parameters,
                              const unsigned int i,
                              libMesh::NumericVector<Real>& sensitivity_rhs);
        
        
        /*!
         *   asks the system to update the nonlinear incompatible mode solution
         */
        void update_incompatible_solution(libMesh::NumericVector<Real>& X,
                                          libMesh::NumericVector<Real>& dX);
 
        
        /*!
         *   Evaluates the volume and boundary outputs for the specified
         *   solution. This reimplements the virtual method from the parent 
         *   class to handle the structural compliance evaluation.
         */
        virtual void calculate_outputs(const libMesh::NumericVector<Real>& X);

        
        /*!
         *   This reimplements the virtual method from the parent
         *   class to handle the structural compliance evaluation.
         */
        void calculate_output_sensitivity(libMesh::ParameterVector& params,
                                          const bool if_total_sensitivity,
                                          const libMesh::NumericVector<Real>& X);

        
    protected:
        
        /*!
         *  calculates the elastic compliance of the system \p S about the
         *  solution defined by \p X. If sensitivity of the quantity is
         *  desired with respect to a parameter, then the parameter can
         *  be specified. If \p dX is provided, then the total sensitivity
         *  is evaluated with respect to the parameter, otherwise the partial
         *  derivative of the output quantity is evaluated.
         */
        virtual void
        _calculate_compliance (const libMesh::NumericVector<Real>& X,
                               libMesh::NonlinearImplicitSystem& S,
                               MAST::RealOutputFunction& output,
                               const MAST::FunctionBase* f = nullptr,
                               const libMesh::NumericVector<Real>* dX = nullptr);
        
        

        /*!
         *   @returns a smart-pointer to a newly created element for
         *   calculation of element quantities.
         */
        virtual std::unique_ptr<MAST::ElementBase>
        _build_elem(const libMesh::Elem& elem);
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void _elem_calculations(MAST::ElementBase& elem,
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
        _elem_linearized_jacobian_solution_product(MAST::ElementBase& elem,
                                                   RealVectorX& vec);

        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                                    bool if_jac,
                                                    RealVectorX& vec,
                                                    RealMatrixX& mat);

        /*!
         *   calculates \f$ d ([J] \{\Delta X\})/ dX  \f$ over \par elem,
         *   and returns the matrix in \par vec .
         */
        virtual void
        _elem_second_derivative_dot_solution_assembly(MAST::ElementBase& elem,
                                                      RealMatrixX& mat);

        /*!
         *   map of local incompatible mode solution per 3D elements
         */
        std::map<const libMesh::Elem*, RealVectorX> _incompatible_sol;
    };
}


#endif // __mast__structural_nonlinear_assembly__
