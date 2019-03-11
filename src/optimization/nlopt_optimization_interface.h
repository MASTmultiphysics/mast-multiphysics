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

#ifndef __MAST_nlopt_optimization_interface_h__
#define __MAST_nlopt_optimization_interface_h__

// MAST includes
#include "optimization/optimization_interface.h"

// NLOpt includes
#include "nlopt.h"

namespace MAST {
    
    class NLOptOptimizationInterface: public MAST::OptimizationInterface {
        
    public:
        
        NLOptOptimizationInterface(nlopt_algorithm alg);
        
        virtual ~NLOptOptimizationInterface()
        { }
        
        virtual void optimize();
        
        /*!
         *   Computes and \returns the value of the objective function for the
         *   specified design variable vector \p x of dimension \p n.
         *   If \p grad is non-null then the gradient of the objective
         *   function is returned in this vector.
         */
        Real
        objective_evaluation(unsigned n,
                             const double* x,
                             double* grad);
        
        
        /*!
         *   Computes the \p m inequality constraints and \returns them in
         *   \p result for the design variable vector \p x of size \p n.
         *   If \p gradient is a non-null vector, then the gradients are
         *   calculated and returned in this vector such that
         *   \f$ \partial f_i/\partial x_j\f$ is stored in \p grad[i*n + j].
         */
        void
        inequality_constraint_evaluation(unsigned m,
                                         double *result,
                                         unsigned n,
                                         const double *x,
                                         double *gradient);
        
    protected:
        
        unsigned int   _iter;
        
        /*!
         *   NLOpt algorithm to use
         */
        nlopt_algorithm _alg;
    };
}





#endif // __MAST_nlopt_optimization_interface_h__
