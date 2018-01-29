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

#ifndef __mast__structural_buckling_eigenproblem_assembly__
#define __mast__structural_buckling_eigenproblem_assembly__


// MAST includes
#include "base/eigenproblem_assembly.h"
#include "base/eigenproblem_assembly_elem_operations.h"

namespace MAST {

    // Forward declerations
    class Parameter;
    
    
    class StructuralBucklingEigenproblemAssembly:
    public MAST::EigenproblemAssembly,
    public MAST::EigenproblemAssemblyElemOperations {
    public:
        
        /*!
         *   constructor associates the eigen system with this assembly object
         */
        StructuralBucklingEigenproblemAssembly();
        
        /*!
         *   destructor resets the association with the eigen system
         *   from this assembly object
         */
        virtual ~StructuralBucklingEigenproblemAssembly();
        
        
        /*!
         *   assembles the global A and B matrices for the modal
         *   eigenvalue problem
         */
        virtual void
        eigenproblem_assemble(libMesh::SparseMatrix<Real>* A,
                              libMesh::SparseMatrix<Real>* B);
        
        
        /*!
         *   set the states and load factors for buckling eigenproblem. If
         *   \p use_linearized_approach is true, then the approximate stability
         *   eigenproblem is solved as \f$ K_1 x = \xi dK x\f$, where
         *   \f$ dK = -(K_2 - K_1) \f$ and \f$ \xi \f$ has the meaning of the
         *   load increment that would lead to a zero eigenvalue of the tangent
         *   stiffness matrix. If the value is false, then the problem is solved
         *   as \f$ K_1 x = \xi (-K_2) x \f$, in which case
         *   \f$ \beta = \frac{\xi}{1+\xi} \f$ is the location between
         *   the load factors \f$\lambda_1 \f$ and \f$\lambda_2 \f$ that 
         *   lead to a zero eigenvalue of the tangent stiffness matrix. In 
         *   both cases the eigenproblem provides an approximation, and the
         *   exact case is satisfied with \f$ \xi = 0 \f$.
         *
         *   Usually it is better to solve the first problem since the 
         *   second case is highly ill-posed and nonlinear around 
         *   \f$ \xi = -1 \f$.
         */
        void set_buckling_data (bool use_linearized_approach,
                                MAST::Parameter& p,
                                const Real lambda1,
                                const Real lambda2,
                                libMesh::NumericVector<Real>& x1,
                                libMesh::NumericVector<Real>& x2);

        
        /*!
         *   calculates the critical load factor based on the eigensolution
         */
        Real critical_point_estimate_from_eigenproblem(Real v) const;
        
        
        /*!
         *   clears the states and load factors for buckling eigenproblem
         */
        virtual void clear_discipline_and_system();
        

        virtual void
        set_elem_solution(const RealVectorX& sol);
        
        /*!
         *   sets the element solution(s) before calculations
         */
        virtual void
        set_elem_sol_sens(MAST::ElementBase& elem,
                          const RealVectorX& sol);
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        elem_calculations(RealMatrixX& mat_A,
                          RealMatrixX& mat_B);
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        elem_sensitivity_calculations(bool base_sol,
                                      RealMatrixX& mat_A,
                                      RealMatrixX& mat_B);
        /*!
         *   initializes the object for the geometric element \p elem. This
         *   expects the object to be in a cleared state, so the user should
         *   call \p clear_elem() between successive initializations.
         */
        virtual void
        init(const libMesh::Elem& elem);
        
        /*!
         *   some simulations frequently deal with 1D/2D elements in 3D space,
         *   which requires use of MAST::LocalElemFE.
         */
        virtual bool
        if_use_local_elem() const {
            
            return true;
        }
        
        /*!
         *   sets additional data for local elem FE.
         */
        virtual void
        set_local_fe_data(MAST::LocalElemFE& fe) const;


    protected:
        
        /*!
         *   map of local incompatible mode solution per 3D elements
         */
        std::map<const libMesh::Elem*, RealVectorX> _incompatible_sol;
        
        
        /*!
         *   whether or not to use the linearized formulation
         */
        bool _use_linearized_formulation;
        
        /*!
         *   load parameter used to define the two stiffness matrices
         */
        MAST::Parameter* _load_param;
        
        
        /*!
         *   values of load factors for which the two stiffness matrices are
         *   calculated for the buckling eigenvalue problem.
         */
        Real _lambda1, _lambda2;
        
        /*!
         *   the equilibrium solutions associated with _lambda1 and _lambda2 
         *   load factors.
         */
        libMesh::NumericVector<Real> *_sol1, *_sol2;
        
        
        /*!
         *   the equilibrium solution sensitivity
         */
        libMesh::NumericVector<Real> *_sol1_sens, *_sol2_sens;
    };
    
}


#endif // __mast__structural_buckling_eigenproblem_assembly__
