/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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

#ifndef __mast__complex_assembly_base_h__
#define __mast__complex_assembly_base_h__


// MAST includes
#include "base/assembly_base.h"

// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"



namespace MAST {
    
    // Forward declerations
    class ComplexSolverBase;
    
    
    class ComplexAssemblyBase:
    public MAST::AssemblyBase,
    public libMesh::NonlinearImplicitSystem::ComputeResidualandJacobian,
    public libMesh::System::SensitivityAssembly {
    public:
        
        /*!
         *   constructor associates this assembly object with the system
         */
        ComplexAssemblyBase();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~ComplexAssemblyBase();

        
        /*!
         *   tells the object to assembly the system of equations for the 
         *   real part when residual_and_jacobian are called.
         */
        void set_assemble_real_part() {
            _if_assemble_real = true;
        }
        
        
        
        /*!
         *   tells the object to assembly the system of equations for the
         *   complex part when residual_and_jacobian are called.
         */
        void set_assemble_imag_part() {
            _if_assemble_real = false;
        }

        
        /*!
         *   calculates the L2 norm of the residual complex system of equations.
         */
        Real residual_l2_norm();
        
        /*!
         *   attaches a system to this discipline, and vice-a-versa
         */
        virtual void
        attach_discipline_and_system(MAST::PhysicsDisciplineBase& discipline,
                                     MAST::SystemInitialization& system)  {
            
            libmesh_error_msg("Error! Invalid function call for ComplexAssemblyBase.");
        }
        
        
        /*!
         *   attaches a system to this discipline, and vice-a-versa
         */
        virtual void
        attach_discipline_and_system(MAST::PhysicsDisciplineBase& discipline,
                                     MAST::ComplexSolverBase& solver,
                                     MAST::SystemInitialization& sys);

        
        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void
        clear_discipline_and_system( );
        
        
        /*!
         *   if the problem is defined about a non-zero base solution,
         *   then this method provides the object with the base solution.
         *   The flag \par if_sens tells the method if \par sol
         *   is the sensitivity of the base solution for the current parameter
         *   being solved for
         */
        void set_base_solution(libMesh::NumericVector<Real>& sol,
                               bool if_sens = false);
        
        /*!
         *   @returns true if a nonzero base solution is used to linearize the
         *   Eigen problem, false otherwise
         */
        bool if_linearized_about_nonzero_solution() const;

        
        /*!
         *   @returns a const reference to the base solution (or
         *   its sensitivity when \par if_sens is true) about which
         *   the Eigen problem was linearized.
         */
        const libMesh::NumericVector<Real>&
        base_sol(bool if_sens = false) const;

        
        /*!
         *    function that assembles the matrices and vectors quantities for
         *    nonlinear solution
         */
        virtual void
        residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                               libMesh::NumericVector<Real>* R,
                               libMesh::SparseMatrix<Real>*  J,
                               libMesh::NonlinearImplicitSystem& S);
        
        void
        residual_and_jacobian_field_split (const libMesh::NumericVector<Real>& X_R,
                                           const libMesh::NumericVector<Real>& X_I,
                                           libMesh::NumericVector<Real>& R_R,
                                           libMesh::NumericVector<Real>& R_I,
                                           libMesh::SparseMatrix<Real>&  J_R,
                                           libMesh::SparseMatrix<Real>&  J_I,
                                           libMesh::NonlinearImplicitSystem& S);

        void
        residual_and_jacobian_blocked (const libMesh::NumericVector<Real>& X,
                                       libMesh::NumericVector<Real>& R,
                                       libMesh::SparseMatrix<Real>&  J,
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
                                        ComplexMatrixX& mat) = 0;
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                                    bool if_jac,
                                                    ComplexVectorX& vec,
                                                    ComplexMatrixX& mat) = 0;
        
                
        /*!
         *   flag tells the solver of the current solution type being sought
         */
        bool _if_assemble_real;


        /*!
         *   pointer to complex solver
         */
        MAST::ComplexSolverBase* _complex_solver;
        
        /*!
         *   base solution about which this problem is defined. This
         *   vector stores the localized values necessary to perform element
         *   calculations.
         */
        libMesh::NumericVector<Real> * _base_sol;
        
        /*!
         *   sensitivity of base solution may be needed for sensitivity
         *   analysis. This vector stores the localized values necessary to
         *   perform element calculations.
         */
        libMesh::NumericVector<Real> * _base_sol_sensitivity;
    };

}



#endif //__mast__complex_assembly_base_h__
