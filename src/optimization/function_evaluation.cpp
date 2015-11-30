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

// MAST includes
#include "optimization/function_evaluation.h"


void
MAST::FunctionEvaluation::output(unsigned int iter, const std::vector<Real> &x,
                                 Real obj, const std::vector<Real> &fval,
                                 bool if_write_to_optim_file) const {
    
    libmesh_assert_equal_to(x.size(), _n_vars);
    libmesh_assert_equal_to(fval.size(), _n_eq + _n_ineq);
    
    
    libMesh::out
    << " *************************** " << std::endl
    << " *** Optimization Output *** " << std::endl
    << " *************************** " << std::endl
    << std::endl
    << "Iter:            "  << std::setw(10) << iter << std::endl
    << "Nvars:           " << std::setw(10) << x.size() << std::endl
    << "Ncons-Equality:  " << std::setw(10) << _n_eq << std::endl
    << "Ncons-Inquality: " << std::setw(10) << _n_ineq << std::endl
    << std::endl
    << "Obj =                  " << std::setw(20) << obj << std::endl
    << std::endl
    << "Vars:            " << std::endl;
    
    for (unsigned int i=0; i<_n_vars; i++)
        libMesh::out
        << "x     [ " << std::setw(10) << i << " ] = "
        << std::setw(20) << x[i] << std::endl;
    
    if (_n_eq) {
        libMesh::out << std::endl
        << "Equality Constraints: " << std::endl;
        
        for (unsigned int i=0; i<_n_eq; i++)
            libMesh::out
            << "feq [ " << std::setw(10) << i << " ] = "
            << std::setw(20) << fval[i] << std::endl;
    }
    
    if (_n_ineq) {
        libMesh::out << std::endl
        << "Inequality Constraints: " << std::endl;
        unsigned int
        n_active      = 0,
        n_violated    = 0,
        max_constr_id = 0;
        Real
        max_constr  = -1.e20;
        
        for (unsigned int i=0; i<_n_ineq; i++) {
            libMesh::out
            << "fineq [ " << std::setw(10) << i << " ] = "
            << std::setw(20) << fval[i+_n_eq];
            if (fabs(fval[i+_n_eq]) <= _tol) {
                n_active++;
                libMesh::out << "  ***" << std::endl;
            }
            else if (fval[i+_n_eq] > _tol) {
                n_violated++;
                libMesh::out << "  +++" << std::endl;
            }
            if (max_constr < fval[i+_n_eq]) {
                max_constr_id = i;
                max_constr    = fval[i+_n_eq];
            }
        }
        
        libMesh::out << std::endl
        << std::setw(35) << " N Active Constraints: "
        << std::setw(20) << n_active << std::endl
        << std::setw(35) << " N Violated Constraints: "
        << std::setw(20) << n_violated << std::endl
        << std::setw(35) << " Most critical constraint: "
        << std::setw(20) << max_constr << std::endl;
    }
    
    libMesh::out << std::endl
    << " *************************** " << std::endl;
    
    
    // the next section writes to the optimization file.
    if (!if_write_to_optim_file)
        return;
    
    {
        // write header for the first iteration
        if (iter == 0) {
            _output << std::setw(10) << "Iter";
            for (unsigned int i=0; i < x.size(); i++) {
                std::stringstream x; x << "x_" << i;
                _output << std::setw(20) << x.str();
            }
            _output << std::setw(20) << "Obj";
            for (unsigned int i=0; i<fval.size(); i++) {
                std::stringstream f; f << "f_" << i;
                _output << std::setw(20) << f.str();
            }
            _output << std::endl;
        }
        
        _output << std::setw(10) << iter;
        for (unsigned int i=0; i < x.size(); i++)
            _output << std::setw(20) << x[i];
        _output << std::setw(20) << obj;
        for (unsigned int i=0; i < fval.size(); i++)
            _output << std::setw(20) << fval[i];
        _output << std::endl;
    }
}





bool
MAST::FunctionEvaluation::verify_gradients(const std::vector<Real>& dvars) {
    
    
    // first call theh evaluate method to get the analytical sensitivities
    Real
    delta           = 1.e-5,
    tol             = 1.e-3,
    obj             = 0.,
    obj_fd          = 0.;
    
    bool
    eval_obj_grad   = true;
    
    
    std::vector<Real>
    dvars_fd   (dvars),
    obj_grad   (_n_vars),
    obj_grad_fd(_n_vars),
    fvals      (_n_ineq + _n_eq),
    fvals_fd   (_n_ineq + _n_eq),
    grads      (_n_vars*(_n_ineq + _n_eq)),
    grads_fd   (_n_vars*(_n_ineq + _n_eq));
    
    
    std::vector<bool>
    eval_grads (_n_eq+_n_ineq);
    
    
    std::fill(    dvars_fd.begin(),     dvars_fd.end(),   0.);
    std::fill(    obj_grad.begin(),     obj_grad.end(),   0.);
    std::fill( obj_grad_fd.begin(),  obj_grad_fd.end(),   0.);
    std::fill(       fvals.begin(),        fvals.end(),   0.);
    std::fill(    fvals_fd.begin(),     fvals_fd.end(),   0.);
    std::fill(       grads.begin(),        grads.end(),   0.);
    std::fill(    grads_fd.begin(),     grads_fd.end(),   0.);
    std::fill(  eval_grads.begin(),   eval_grads.end(), true);
    
    
    // calculate the analytical sensitivity
    this->evaluate(dvars,
                   obj,
                   eval_obj_grad,
                   obj_grad,
                   fvals,
                   eval_grads,
                   grads);
    
    
    // now turn off the sensitivity variables
    eval_obj_grad = false;
    std::fill(  eval_grads.begin(),   eval_grads.end(), false);
    
    // now iteratve over the design variables, and calculate the finite
    // difference sensitivity values
    
    for (unsigned int i=0; i<_n_vars; i++) {
        
        // copy the original vector
        dvars_fd =  dvars;
        
        // now perturb it
        dvars_fd[i] += delta;
        
        // call the evaluate routine
        obj_fd       = 0.;
        std::fill(    fvals_fd.begin(),     fvals_fd.end(),   0.);
        this->evaluate(dvars_fd,
                       obj_fd,
                       eval_obj_grad,
                       obj_grad_fd,
                       fvals_fd,
                       eval_grads,
                       grads_fd);
        
        // objective gradient
        obj_grad_fd[i]  = (obj_fd-obj)/delta;
        
        // constraint gradient
        for (unsigned int j=0; j<_n_eq+_n_ineq; j++)
            grads_fd[i*(_n_eq+_n_ineq)+j]  = (fvals_fd[j]-fvals[j])/delta;
    }
    
    
    // compare the values
    std::cout
    << " *** Objective function gradients: analytical vs numerical"
    << std::endl;
    
    bool accurate_sens = true;
    
    for (unsigned int i=0; i<_n_vars; i++)
        if (fabs(obj_grad[i] - obj_grad_fd[i])/obj_grad[i] > tol) {
            std::cout
            << " Mismatched sensitivity: DV:  "  << i << "   "
            << obj_grad[i] << "    " << obj_grad_fd[i] << std::endl;
            accurate_sens = false;
        }
    
    
    
    std::cout
    << " *** Constraint function gradients: analytical vs numerical"
    << std::endl;
    
    for (unsigned int j=0; j<_n_eq+_n_ineq; j++) {
        std::cout << "  Constraint: " << j << std::endl;
        for (unsigned int i=0; i<_n_vars; i++)
            if (fabs(grads[i*(_n_eq+_n_ineq)+j] - grads_fd[i*(_n_eq+_n_ineq)+j])/grads[i*(_n_eq+_n_ineq)+j] > tol) {
                
                std::cout
                << " Mismatched sensitivity:  DV:  "  << i << "   "
                << grads[i*(_n_eq+_n_ineq)+j] << "    "
                << grads_fd[i*(_n_eq+_n_ineq)+j] << std::endl;
                accurate_sens = false;
            }
    }
    
    // print the message that all sensitivity data satisfied limits.
    if (accurate_sens)
        std::cout
        << "Verify gradients: all gradients satisfied relative tol: " << tol
        << "  with delta:  " << delta
        << std::endl;
    
    return accurate_sens;
}



