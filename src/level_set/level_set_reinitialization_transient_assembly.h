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

#ifndef __mast__level_set_reinitialization_transient_assembly__
#define __mast__level_set_reinitialization_transient_assembly__

// MAST includes
#include "base/transient_assembly.h"


namespace MAST {

    class LevelSetReinitializationTransientAssembly:
    public MAST::TransientAssembly {
      
    public:
        
        LevelSetReinitializationTransientAssembly();
        
        virtual ~LevelSetReinitializationTransientAssembly();

        /*!
         *   For reinitialization to \f$ |\nabla(\phi)| = 1 \f$, the solution
         *   before initialization is used to calculate the source and velocity
         *   switching. This method sets that solution.
         */
        void set_reference_solution(const libMesh::NumericVector<Real>& sol);

        void clear_reference_solution(const libMesh::NumericVector<Real>& sol) {
            
            _ref_sol = nullptr;
        }

        /*!
         *    function that assembles the matrices and vectors quantities for
         *    nonlinear solution
         */
        virtual void
        residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                               libMesh::NumericVector<Real>* R,
                               libMesh::SparseMatrix<Real>*  J,
                               libMesh::NonlinearImplicitSystem& S);

        
    protected:
        
        const libMesh::NumericVector<Real>* _ref_sol;
    };
}


#endif // __mast__level_set_reinitialization_transient_assembly__
