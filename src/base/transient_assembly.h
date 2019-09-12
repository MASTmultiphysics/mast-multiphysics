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

#ifndef __mast__transient_assembly__
#define __mast__transient_assembly__

// MAST includes
#include "base/assembly_base.h"

// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"



namespace MAST {
    
    // Forward declerations
    class TransientSolverBase;
    class TransientAssemblyElemOperations;
    
    
    class TransientAssembly:
    public MAST::AssemblyBase {
    public:
        
        /*!
         *    user-provided object to perform actions
         *    after assembly and before returning to the solver. Use
         *    \p set_post_assembly_object to provide a pointer to the object.
         */
        class PostAssemblyOperation{
            
        public:
            PostAssemblyOperation() {}
            virtual ~PostAssemblyOperation() {}
            virtual void post_assembly(const libMesh::NumericVector<Real>& X,
                                       libMesh::NumericVector<Real>* R,
                                       libMesh::SparseMatrix<Real>*  J,
                                       libMesh::NonlinearImplicitSystem& S) = 0;
        };

        /*!
         *   constructor associates this assembly object with the system
         */
        TransientAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~TransientAssembly();
        
        
        /*!
         *    sets the PostAssemblyOperation object for use after assembly.
         *    Note that calling \p clear_discipline_and_system() will
         *    clear this pointer and the user will have to call this function
         *    again.
         */
        void
        set_post_assembly_operation(MAST::TransientAssembly::PostAssemblyOperation& post);


        /*!
         *    function that assembles the matrices and vectors quantities for
         *    nonlinear solution
         */
        virtual void
        residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                               libMesh::NumericVector<Real>* R,
                               libMesh::SparseMatrix<Real>*  J,
                               libMesh::NonlinearImplicitSystem& S);
        
        /*!
         *    calculates the product of the Jacobian and a perturbation in solution
         *    vector \f$ [J] \{\Delta X\}  \f$. For a single discipline system the
         *    solution vector and linearized solution provided here are used. For
         *    a multiphysics system, the user must ensure that all relevant
         *    multidisciplinary data-structures are initialized before calling
         *    this method.
         */
        virtual void
        linearized_jacobian_solution_product(const libMesh::NumericVector<Real>& X,
                                             const libMesh::NumericVector<Real>& dX,
                                             libMesh::NumericVector<Real>& JdX,
                                             libMesh::NonlinearImplicitSystem& S);
        
        /**
         * Assembly function.  This function will be called
         * to assemble the sensitivity of system residual prior to a solve and must
         * be provided by the user in a derived class. The method provides dR/dp_i
         * for \p i ^th parameter in the vector \p parameters.
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
         *    this object, if non-NULL is user-provided to perform actions
         *    after assembly and before returning to the solver
         */
        MAST::TransientAssembly::PostAssemblyOperation* _post_assembly;

    };
    
    
}

#endif // __mast__transient_assembly__
