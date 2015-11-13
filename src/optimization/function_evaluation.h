/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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


#ifndef __mast__function_evaluation_h__
#define __mast__function_evaluation_h__

// C++ includes
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>


// MAST includes
#include "base/mast_data_types.h"


namespace MAST {
    
    class FunctionEvaluation {
        
    public:
        
        FunctionEvaluation(std::ostream& tab_output):
        _n_vars(0),
        _n_eq(0),
        _n_ineq(0),
        _max_iters(0),
        _n_rel_change_iters(5),
        _tol(1.0e-6),
        _output(tab_output)
        { }
        
        virtual ~FunctionEvaluation() { }
        
        
        unsigned int n_vars() const {
            return _n_vars;
        }
        
        
        unsigned int n_eq() const {
            return _n_eq;
        }
        
        
        unsigned int n_ineq() const{
            return _n_ineq;
        }
        
        
        unsigned int max_iters() const{
            return _max_iters;
        }
        
        
        unsigned int n_iters_relative_change() const {
            return _n_rel_change_iters;
        }
        
        
        Real tolerance() const{
            return _tol;
        }
        
        
        virtual void init_dvar(std::vector<Real>& x,
                               std::vector<Real>& xmin,
                               std::vector<Real>& xmax) = 0;
        
        /*!
         *   \par grads(k): Derivative of f_i(x) with respect
         *   to x_j, where k = (j-1)*M + i.
         */
        virtual void evaluate(const std::vector<Real>& dvars,
                              Real& obj,
                              bool eval_obj_grad,
                              std::vector<Real>& obj_grad,
                              std::vector<Real>& fvals,
                              std::vector<bool>& eval_grads,
                              std::vector<Real>& grads) = 0;
        
        virtual void output(unsigned int iter,
                            const std::vector<Real>& x,
                            Real obj,
                            const std::vector<Real>& fval,
                            bool if_write_to_optim_file) const;
        
        
        /*!
         *  verifies the gradients at the specified design point
         */
        virtual bool verify_gradients(const std::vector<Real>& dvars);
        
        
        
    protected:
        
        unsigned int _n_vars;
        
        unsigned int _n_eq;
        
        unsigned int _n_ineq;
        
        unsigned int _max_iters;
        
        unsigned int _n_rel_change_iters;
        
        Real _tol;
        
        std::ostream& _output;
    };


}

#endif // __mast__function_evaluation_h__

