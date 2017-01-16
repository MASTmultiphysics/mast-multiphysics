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

#ifndef __mast__linearized_assembly_base__
#define __mast__linearized_assembly_base__

// MAST includes
#include "base/nonlinear_implicit_assembly.h"



namespace MAST {

    // Forward declerations
    class Parameter;

    
    /*!
     *   Assembles the system of equations for a problem linearized about a 
     *   base solution. The base solution is provided to this assembly class
     *   and the solution vector provided in the assemble function is the 
     *   perturbation about this solution. While the \p EigenproblemAssembly and
     *   \p ComplexAssembly classes also provide a similar linearized assembly, 
     *   the focus in the former class is for assembly of coefficient matrices,
     *   and the assembly of complex residual and Jacobian quantities in the 
     *   latter class. In this class, the focus is on calculation of real-valued
     *   time-domain residual and Jacobian quantities. While the Jacobian can be
     *   same as that provided by the NonlinearImplicitAssembly class due to its 
     *   dependence only on the base solution (except in cases where additional
     *   stabilization terms are added to the linearized system, eg. to
     *   handle discontinuities in supersonic linearized flow), the residual
     *   accounts for both the base and perturbed solution.
     */
    
    class LinearizedAssemblyBase:
    public MAST::NonlinearImplicitAssembly {
        
    public:
        
        /*!
         *   constructor associates the eigen system with this assembly object
         */
        LinearizedAssemblyBase();
        
        /*!
         *   destructor resets the association with the eigen system
         *   from this assembly object
         */
        virtual ~LinearizedAssemblyBase();
        
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
         *   if the eigenproblem is defined about a non-zero base solution,
         *   then this method provides the object with the base solution.
         *   The flag \par if_sens tells the method if \par sol
         *   is the sensitivity of the base solution for the current parameter
         *   being solved for
         */
        void set_base_solution(const libMesh::NumericVector<Real>& sol,
                               bool if_sens = false);
        
        
        /*!
         *   Clears the pointer to the solution.
         *   The flag \par if_sens tells the method to clear the pointer
         *   the solution sensitivity vector.
         */
        void clear_base_solution(bool if_sens = false);
        
        
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
         *   @returns a non-const reference to the base solution (or
         *   its sensitivity when \par if_sens is true) about which
         *   the Eigen problem was linearized.
         */
        libMesh::NumericVector<Real>& base_sol(bool if_sens = false);
        
        
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
        
    protected:
        
        
        /*!
         *   base solution about which this eigenproblem is defined. This
         *   vector stores the localized values necessary to perform element
         *   calculations.
         */
        const libMesh::NumericVector<Real> * _base_sol;
        
        /*!
         *   sensitivity of base solution may be needed for sensitivity
         *   analysis. This vector stores the localized values necessary to
         *   perform element calculations.
         */
        const libMesh::NumericVector<Real> * _base_sol_sensitivity;
        
    };
    
}

#endif // __mast__linearized_assembly_base__
