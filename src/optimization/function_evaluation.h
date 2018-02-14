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


#ifndef __mast__function_evaluation_h__
#define __mast__function_evaluation_h__

// C++ includes
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>


// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/parallel_object.h"


namespace MAST {
    
    class FunctionEvaluation:
    public libMesh::ParallelObject {
        
    public:
        
        FunctionEvaluation(const libMesh::Parallel::Communicator& comm_in):
        libMesh::ParallelObject(comm_in),
        _n_vars(0),
        _n_eq(0),
        _n_ineq(0),
        _max_iters(0),
        _n_rel_change_iters(5),
        _tol(1.0e-6),
        _output(nullptr)
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
        
        
        /*!
         *   sets the output file and the function evaluation will 
         *   write the optimization iterates to this file. If this is not called
         *   no file will be created. The user may choose to set different 
         *   file names on different ranks on the communicator, or set the 
         *   output only on one of the ranks. Specifying the same filename 
         *   on multiple ranks will lead to undefined behavior since all ranks
         *   will try to write to the same file.
         */
        void set_output_file(const std::string& nm) {
            
            if (!_output)
                _output = new std::ofstream;
            
            _output->close();
            _output->open(nm.c_str(), std::ofstream::out);
        }

        
        /*!
         *   outputs the the current iterate to libMesh::out, and to the 
         *   output file if it was set for this rank.
         */
        virtual void output(unsigned int iter,
                            const std::vector<Real>& x,
                            Real obj,
                            const std::vector<Real>& fval,
                            bool if_write_to_optim_file) const;
        
        
        /*!
         *   This reads and initializes the DV vector from a previous
         *   optimization history output file. This will verify that the
         *   optimization setup (number of DVs, constraints, etc.) are the
         *   same as the initialized data for this object and then
         *   read the dv values from \par iter iteration into
         *   \par x. 
         */
        void initialize_dv_from_output_file(const std::string& nm,
                                            const unsigned int iter,
                                            std::vector<Real> &x);
        
        /*!
         *  verifies the gradients at the specified design point
         */
        virtual bool verify_gradients(const std::vector<Real>& dvars);
        
        
#if MAST_ENABLE_NPSOL == 1
        typedef void (*funobj) (int*    mode,
                                int*    n,
                                double* x,
                                double* f,
                                double* g,
                                int*    nstate);
        
        
        typedef void (*funcon) (int*    mode,
                                int*    ncnln,
                                int*    n,
                                int*    ldJ,
                                int*    needc,
                                double* x,
                                double* c,
                                double* cJac,
                                int*    nstate);

        /*!
         *  @returns a pointer to the function that evaluates the objective
         *  used for NPSOL interface
         */
        virtual funobj
        get_objective_evaluation_function() {
            
            // should not get here, if the derived method implements its
            // specialized method
            libmesh_assert(false);
            return nullptr;
        }

        
        /*!
         *  @returns a pointer to the function that evaluates the constraint
         *  used for NPSOL interface
         */
        virtual funcon
        get_constraint_evaluation_function() {
            
            // should not get here, if the derived method implements its
            // specialized method
            libmesh_assert(false);
            return nullptr;
        }
#endif
        
        
    protected:
        
        unsigned int _n_vars;
        
        unsigned int _n_eq;
        
        unsigned int _n_ineq;
        
        unsigned int _max_iters;
        
        unsigned int _n_rel_change_iters;
        
        Real _tol;
        
        std::ofstream* _output;
    };


}

#endif // __mast__function_evaluation_h__

