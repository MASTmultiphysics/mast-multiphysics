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

#ifndef __mast__generalized_alpha_transient_solver__
#define __mast__generalized_alpha_transient_solver__

// MAST includes
#include "solver/second_order_newmark_transient_solver.h"


namespace MAST {
    
    
    /*!
     *    This class implements the generalized alpha method for solution of a
     *    second-order ODE.
     *
     */
    class GeneralizedAlphaTransientSolver:
    public MAST::SecondOrderNewmarkTransientSolver {
    public:
        GeneralizedAlphaTransientSolver();
        
        virtual ~GeneralizedAlphaTransientSolver();
        
        /*!
         *    \f$ \rho_\infty \f$ parameter used by this solver. A value of 1 provides no
         *    dissippation and a value of 0 provides asymptotic annihilation. Default is 0.2.
         */
        Real rho_inf;

        /*!
         *    \f$ \alpha_m \f$ parameter used by this solver.
         */
        Real alpha_m;

        /*!
         *    \f$ \alpha_f \f$ parameter used by this solver.
         */
        Real alpha_f;

        /*!
         *    Computes the value of coefficients basde on value of \f$ \rho_\infty \f$
         */
        void update_coefficient(Real rho_infinity);

    
        /*!
         *    provides the element with the transient data for calculations
         */
        virtual void
        set_element_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                         const std::vector<libMesh::NumericVector<Real>*>& sols);
        
        /*!
         *    provides the element with the sensitivity of transient data for
         *    calculations
         */
        virtual void
        extract_element_sensitivity_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                                         const std::vector<libMesh::NumericVector<Real>*>& sols,
                                         std::vector<RealVectorX>& local_sols);

        /*!
         *    provides the element with the transient data for calculations
         */
        virtual void
        set_element_perturbed_data
        (const std::vector<libMesh::dof_id_type>& dof_indices,
         const std::vector<libMesh::NumericVector<Real>*>& sols);
        
        
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
        
    protected:
        
    };
    
}

#endif // __mast__generalized_alpha_transient_solver__
