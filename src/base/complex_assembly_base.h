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

#ifndef __mast__complex_assembly_base_h__
#define __mast__complex_assembly_base_h__


// MAST includes
#include "base/assembly_base.h"

// libMesh includes
#include "libmesh/libmesh_logging.h"
#include "libmesh/nonlinear_implicit_system.h"



namespace MAST {
    
    // Forward declerations
    class ComplexSolverBase;
    class Parameter;
    
    
    class ComplexAssemblyBase:
    public MAST::AssemblyBase {
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
         *   calculates the L2 norm of the residual complex system of equations.
         */
        Real residual_l2_norm(const libMesh::NumericVector<Real>& real,
                              const libMesh::NumericVector<Real>& imag);
        
        
        /*!
         *   if the problem is defined about a non-zero base solution,
         *   then this method provides the object with the base solution.
         *   The flag \p if_sens tells the method if \p sol
         *   is the sensitivity of the base solution for the current parameter
         *   being solved for
         */
        void set_base_solution(const libMesh::NumericVector<Real>& sol,
                               bool if_sens = false);

        
        /*!
         *   Clears pointer to the solution vector
         *   The flag \p if_sens tells the method to clear pointer of the 
         *   solution sensitivity
         */
        void clear_base_solution(bool if_sens = false);

        
        /*!
         *   @returns true if a nonzero base solution is used to linearize the
         *   Eigen problem, false otherwise
         */
        bool if_linearized_about_nonzero_solution() const;

        
        /*!
         *   @returns a const reference to the base solution (or
         *   its sensitivity when \p if_sens is true) about which
         *   the Eigen problem was linearized.
         */
        const libMesh::NumericVector<Real>& base_sol(bool if_sens = false) const;
        
        
        void
        residual_and_jacobian_field_split (const libMesh::NumericVector<Real>& X_R,
                                           const libMesh::NumericVector<Real>& X_I,
                                           libMesh::NumericVector<Real>& R_R,
                                           libMesh::NumericVector<Real>& R_I,
                                           libMesh::SparseMatrix<Real>&  J_R,
                                           libMesh::SparseMatrix<Real>&  J_I);

        /*!
         *   Assembles the residual and Jacobian of the N complex system of equations
         *   split into 2N real system of equations. The solution vector \p X
         *   is the current complex solution with the real and imaginary parts 
         *   of each element stored as adjacent entries. Likewise, the Jaacobian
         *   matrix has a 2x2 block storage. If \p p is provided, then \p R
         *   will return the sensitivity of the residual vector.
         */
        void
        residual_and_jacobian_blocked (const libMesh::NumericVector<Real>& X,
                                       libMesh::NumericVector<Real>& R,
                                       libMesh::SparseMatrix<Real>&  J,
                                       MAST::Parameter* p = nullptr);

        /**
         * Assembly function.  This function will be called
         * to assemble the RHS of the sensitivity equations (which is -1 times
         * sensitivity of system residual) prior to a solve and must
         * be provided by the user in a derived class. The method provides dR/dp
         * for \p f parameter.
         *
         * If the routine is not able to provide sensitivity for this parameter,
         * then it should return false, and the system will attempt to use
         * finite differencing.
         */
        virtual bool
        sensitivity_assemble (const MAST::FunctionBase& f,
                              libMesh::NumericVector<Real>& sensitivity_rhs);
        
    protected:
        
        /*!
         *   base solution about which this problem is defined. This
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



#endif //__mast__complex_assembly_base_h__
