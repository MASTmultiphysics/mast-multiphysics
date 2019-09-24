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

// C++ includes
#include <sys/stat.h>
#include <string>
#include <boost/algorithm/string.hpp>

// MAST includes
#include "optimization/function_evaluation.h"

// libMesh includes
#include "libmesh/parallel_implementation.h"


void
MAST::FunctionEvaluation::attach_optimization_interface(MAST::OptimizationInterface& opt) {
    
    libmesh_assert(!_optimization_interface);
    
    _optimization_interface = &opt;
}


void
MAST::FunctionEvaluation::output(unsigned int iter,
                                 const std::vector<Real> &x,
                                 Real obj,
                                 const std::vector<Real> &fval,
                                 bool if_write_to_optim_file) {
    
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
                libMesh::out << "  ***";
            }
            else if (fval[i+_n_eq] > _tol) {
                n_violated++;
                libMesh::out << "  +++";
            }
            libMesh::out  << std::endl;
            
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
    if (!if_write_to_optim_file ||
        !_output)  // or if the output has not been specified
        return;
    
    {
        // write header for the first iteration
        if (iter == 0) {

            // number of desing variables
            *_output
            << std::setw(10) << "n_dv" << std::setw(10) << _n_vars << std::endl;
            *_output
            << std::setw(10) << "n_eq" << std::setw(10) << _n_eq << std::endl;
            *_output
            << std::setw(10) << "n_ineq" << std::setw(10) << _n_ineq << std::endl;

            *_output << std::setw(10) << "Iter";
            for (unsigned int i=0; i < x.size(); i++) {
                std::stringstream x; x << "x_" << i;
                *_output << std::setw(20) << x.str();
            }
            *_output << std::setw(20) << "Obj";
            for (unsigned int i=0; i<fval.size(); i++) {
                std::stringstream f; f << "f_" << i;
                *_output << std::setw(20) << f.str();
            }
            *_output << std::endl;
        }
        
        *_output << std::setw(10) << iter;
        for (unsigned int i=0; i < x.size(); i++)
            *_output << std::setw(20) << x[i];
        *_output << std::setw(20) << obj;
        for (unsigned int i=0; i < fval.size(); i++)
            *_output << std::setw(20) << fval[i];
        *_output << std::endl;
    }
}



void
MAST::FunctionEvaluation::initialize_dv_from_output_file(const std::string& nm,
                                                         const unsigned int iter,
                                                         std::vector<Real> &x) {
    
    struct stat stat_info;
    int stat_result = stat(nm.c_str(), &stat_info);
    
    if (stat_result != 0)
        libmesh_error_msg("File does not exist: " + nm);
    
    if (!std::ifstream(nm))
        libmesh_error_msg("File missing: " + nm);
    
    std::ifstream input;
    input.open(nm, std::ofstream::in);
    
    
    std::string
    line;
    unsigned int
    ndv        = 0,
    nineq      = 0,
    neq        = 0,
    it_num     = 0;
    
    std::vector<std::string> results;
    
    // number of desing variables
    std::getline(input, line);
    boost::trim(line);
    boost::split(results, line, boost::is_any_of(" \t"), boost::token_compress_on);
    libmesh_assert_equal_to(results[0],   "n_dv");
    ndv = stod(results[1]);
    libmesh_assert_equal_to(  ndv, x.size());
    
    
    // number of equality constraint
    std::getline(input, line);
    boost::trim(line);
    boost::split(results, line, boost::is_any_of(" \t"), boost::token_compress_on);
    libmesh_assert_equal_to(results[0],   "n_eq");
    neq = stod(results[1]);
    libmesh_assert_equal_to(  neq, _n_eq);
    
    
    // number of inequality constriants
    std::getline(input, line);
    boost::trim(line);
    boost::split(results, line, boost::is_any_of(" \t"), boost::token_compress_on);
    libmesh_assert_equal_to(results[0],   "n_ineq");
    nineq = stod(results[1]);
    //libmesh_assert_equal_to(  nineq, _n_ineq);
    
    
    // skip all lines before iter.
    while (!input.eof() && it_num < iter+1) {
        std::getline(input, line);
        it_num++;
    }
    
    // make sure that the iteration number is what we are looking for
    std::getline(input, line);
    boost::trim(line);
    boost::split(results, line, boost::is_any_of(" \t"), boost::token_compress_on);
    
    libmesh_assert_greater(results.size(), ndv+1);
    
    it_num = stoi(results[0]);
    libmesh_assert_equal_to(it_num, iter);
    
    // make sure that the file has data
    for (unsigned int i=0; i<ndv; i++)
        x[i] = stod(results[i+1]);
}


bool
MAST::FunctionEvaluation::verify_gradients(const std::vector<Real>& dvars) {
    
    
    // first call theh evaluate method to get the analytical sensitivities
    Real
    delta           = 1.e-5,
    tol             = 1.e-3,
    obj             = 0.,
    obj_fd_p        = 0.,  // at x+h
    obj_fd_m        = 0.;  // at x-h

    bool
    eval_obj_grad   = true;
    
    
    std::vector<Real>
    dvars_fd   (dvars),
    obj_grad   (_n_vars),
    obj_grad_fd(_n_vars),
    fvals      (_n_ineq + _n_eq),
    fvals_fd_p (_n_ineq + _n_eq),  // at x+h
    fvals_fd_m (_n_ineq + _n_eq),  // at x-h
    grads      (_n_vars*(_n_ineq + _n_eq)),
    grads_fd   (_n_vars*(_n_ineq + _n_eq));
    
    
    std::vector<bool>
    eval_grads (_n_eq+_n_ineq);
    
    
    std::fill(    dvars_fd.begin(),     dvars_fd.end(),   0.);
    std::fill(    obj_grad.begin(),     obj_grad.end(),   0.);
    std::fill( obj_grad_fd.begin(),  obj_grad_fd.end(),   0.);
    std::fill(       fvals.begin(),        fvals.end(),   0.);
    std::fill(  fvals_fd_p.begin(),   fvals_fd_p.end(),   0.);
    std::fill(  fvals_fd_m.begin(),   fvals_fd_m.end(),   0.);
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
        
        // central difference approx
        // du/dx = ((u2-u1)/(x2-x1) + (u1-u0)/(x1-x0))/2
        //       = ((u2-u1)/h + (u1-u0)/h)/2
        //       = (u2-u0)/2h
        
        // now perturb it: first the positive value
        dvars_fd[i] += delta;
        
        // call the evaluate routine
        obj_fd_p     = 0.;
        obj_fd_m     = 0.;
        std::fill(  fvals_fd_p.begin(),   fvals_fd_p.end(),   0.);
        std::fill(  fvals_fd_m.begin(),   fvals_fd_m.end(),   0.);
        this->evaluate(dvars_fd,
                       obj_fd_p,
                       eval_obj_grad,
                       obj_grad_fd,
                       fvals_fd_p,
                       eval_grads,
                       grads_fd);

        // now perturb it: first the positive value
        dvars_fd[i] -= 2*delta;
        
        this->evaluate(dvars_fd,
                       obj_fd_m,
                       eval_obj_grad,
                       obj_grad_fd,
                       fvals_fd_m,
                       eval_grads,
                       grads_fd);

        // objective gradient
        obj_grad_fd[i]  = (obj_fd_p-obj_fd_m)/2./delta;
        
        // constraint gradient
        for (unsigned int j=0; j<_n_eq+_n_ineq; j++)
            grads_fd[i*(_n_eq+_n_ineq)+j]  = (fvals_fd_p[j]-fvals_fd_m[j])/2./delta;
    }
    
    
    // compare the values
    libMesh::out
    << " *** Objective function gradients: analytical vs numerical"
    << std::endl;
    
    bool accurate_sens = true;

    libMesh::out
    << std::setw(10) << "DV"
    << std::setw(30) << "Analytical"
    << std::setw(30) << "Numerical" << std::endl;

    for (unsigned int i=0; i<_n_vars; i++) {
        libMesh::out
        << std::setw(10) << i
        << std::setw(30) << obj_grad[i]
        << std::setw(30) << obj_grad_fd[i];
        if (fabs((obj_grad[i] - obj_grad_fd[i])/obj_grad[i]) > tol) {
            libMesh::out << " : Mismatched sensitivity";
            accurate_sens = false;
        }
        libMesh::out << std::endl;
    }
    
    
    
    libMesh::out
    << " *** Constraint function gradients: analytical vs numerical"
    << std::endl;

    libMesh::out
    << std::setw(10) << "DV"
    << std::setw(30) << "Analytical"
    << std::setw(30) << "Numerical" << std::endl;

    for (unsigned int j=0; j<_n_eq+_n_ineq; j++) {
        
        libMesh::out << "  Constraint: " << j << std::endl;
        for (unsigned int i=0; i<_n_vars; i++) {
            libMesh::out
            << std::setw(10) << i
            << std::setw(30) << grads[i*(_n_eq+_n_ineq)+j]
            << std::setw(30) << grads_fd[i*(_n_eq+_n_ineq)+j];
            if (fabs((grads[i*(_n_eq+_n_ineq)+j] - grads_fd[i*(_n_eq+_n_ineq)+j])/grads[i*(_n_eq+_n_ineq)+j]) > tol) {
                libMesh::out << " : Mismatched sensitivity";
                accurate_sens = false;
            }
            libMesh::out << std::endl;
        }
    }
    // print the message that all sensitivity data satisfied limits.
    if (accurate_sens)
        libMesh::out
        << "Verify gradients: all gradients satisfied relative tol: " << tol
        << "  with delta:  " << delta
        << std::endl;
    
    return accurate_sens;
}


void
MAST::FunctionEvaluation::parametric_line_study(const std::string& nm,
                                                const unsigned int iter1,
                                                const unsigned int iter2,
                                                unsigned int divs) {

    std::vector<Real>
    dv1(_n_vars, 0.),
    dv2(_n_vars, 0.),
    dv (_n_vars, 0.),
    fval(_n_ineq+_n_eq, 0.);
    
    this->initialize_dv_from_output_file(nm, iter1, dv1);
    this->initialize_dv_from_output_file(nm, iter2, dv2);
    
    Real
    f   = 0.,
    obj = 0.;

    for (unsigned int i=0; i<=divs; i++) {
        
        f = (1.*i)/(1.*divs);
        for (unsigned int j=0; j<_n_vars; j++) {
            dv[j] = (1.-f) * dv1[j] + f * dv2[j];
        }
        
        this->_output_wrapper(i, dv, obj, fval, true);
    }
}



void
MAST::FunctionEvaluation::sanitize_parallel() {

    unsigned int
    N                  = this->n_vars(),
    N_EQ               = this->n_eq(),
    N_INEQ             = this->n_ineq(),
    n_rel_change_iters = this->n_iters_relative_change();
    
    // make sure all processors have the same values
    libmesh_assert(this->comm().verify(N));
    libmesh_assert(this->comm().verify(N_EQ));
    libmesh_assert(this->comm().verify(N_INEQ));
    libmesh_assert(this->comm().verify(n_rel_change_iters));
}



void
MAST::FunctionEvaluation::_init_dvar_wrapper(std::vector<Real>& x,
                                             std::vector<Real>& xmin,
                                             std::vector<Real>& xmax) {

    this->init_dvar(x, xmin, xmax);
    
    libmesh_assert(this->comm().verify(x));
    libmesh_assert(this->comm().verify(xmin));
    libmesh_assert(this->comm().verify(xmax));
}



void
MAST::FunctionEvaluation::_evaluate_wrapper(const std::vector<Real>& dvars,
                                            Real& obj,
                                            bool eval_obj_grad,
                                            std::vector<Real>& obj_grad,
                                            std::vector<Real>& fvals,
                                            std::vector<bool>& eval_grads,
                                            std::vector<Real>& grads) {
    
    // verify that all values going into the function are consistent
    // across all processors
    libmesh_assert(this->comm().verify(dvars));
    libmesh_assert(this->comm().verify(eval_obj_grad));
    
    this->evaluate(dvars,
                   obj,
                   eval_obj_grad,
                   obj_grad,
                   fvals,
                   eval_grads,
                   grads);
    
    // verify that all output values coming out of all functions are
    // consistent across all processors
    libmesh_assert(this->comm().verify(obj));
    libmesh_assert(this->comm().verify(obj_grad));
    libmesh_assert(this->comm().verify(fvals));
    libmesh_assert(this->comm().verify(grads));
}



void
MAST::FunctionEvaluation::_output_wrapper(unsigned int iter,
                                          const std::vector<Real>& x,
                                          Real obj,
                                          const std::vector<Real>& fval,
                                          bool if_write_to_optim_file) {
    
    // verify that all values going into the function are consistent
    // across all processors
    libmesh_assert(this->comm().verify(iter));
    libmesh_assert(this->comm().verify(x));
    
    this->output(iter, x, obj, fval, if_write_to_optim_file);
}

