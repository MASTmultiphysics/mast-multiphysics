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

#ifndef __mast__structural_buckling_eigenproblem_assembly__
#define __mast__structural_buckling_eigenproblem_assembly__


// MAST includes
#include "base/eigenproblem_assembly.h"


namespace MAST {

    // Forward declerations
    class Parameter;
    
    
    class StructuralBucklingEigenproblemAssembly:
    public MAST::EigenproblemAssembly {
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
        
        /**
         * Assembly function.  This function will be called
         * to assemble the sensitivity of eigenproblem matrices.
         * The method provides dA/dp_i and dB/dpi for \par i ^th parameter
         * in the vector \par parameters.
         *
         * If the routine is not able to provide sensitivity for this parameter,
         * then it should return false, and the system will attempt to use
         * finite differencing.
         */
        virtual bool
        eigenproblem_sensitivity_assemble (const libMesh::ParameterVector& parameters,
                                           const unsigned int i,
                                           libMesh::SparseMatrix<Real>* sensitivity_A,
                                           libMesh::SparseMatrix<Real>* sensitivity_B);
        
        
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
        
        
        /*!
         *   @returns a MAST::FEBase object for calculation of finite element
         *   quantities. This creates LocalElemFE for 1D and 2D elements.
         */
        virtual std::unique_ptr<MAST::FEBase>
        build_fe(const libMesh::Elem& e);


    protected:
        
        /*!
         *   @returns a smart-pointer to a newly created element for
         *   calculation of element quantities.
         */
        virtual std::unique_ptr<MAST::ElementBase>
        _build_elem(const libMesh::Elem& elem);
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        _elem_calculations(MAST::ElementBase& elem,
                           RealMatrixX& mat_A,
                           RealMatrixX& mat_B) {
            libmesh_error(); // should not get here.
        }
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                       RealMatrixX& mat_A,
                                       RealMatrixX& mat_B) {
            libmesh_error(); // should not get here.
        }

        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        _elem_calculations(MAST::ElementBase& elem,
                           RealMatrixX& mat_A);
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element matrices for the eigenproblem
         *   \f$ A x = \lambda B x \f$.
         */
        virtual void
        _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                       RealMatrixX& mat_A);
        
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
