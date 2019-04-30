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

#ifndef __MAST_npsol_optimization_interface_h__
#define __MAST_npsol_optimization_interface_h__

// C++ includes
#include <map>

// MAST includes
#include "optimization/optimization_interface.h"



namespace MAST {
    
    class NPSOLOptimizationInterface: public MAST::OptimizationInterface {
        
    public:
        
        NPSOLOptimizationInterface();
        
        virtual ~NPSOLOptimizationInterface() { }
        
        virtual void optimize();
        
        virtual void
        attach_function_evaluation_object (MAST::FunctionEvaluation& feval);

    protected:

        void  _print_termination_message(const int INFORM);

        
        void (*_funobj) (int*    mode,
                         int*    n,
                         double* x,
                         double* f,
                         double* g,
                         int*    nstate);
        
        void (*_funcon) (int*    mode,
                         int*    ncnln,
                         int*    n,
                         int*    ldJ,
                         int*    needc,
                         double* x,
                         double* c,
                         double* cJac,
                         int*    nstate);

        std::map<int, std::string> _exit_message;
        std::map<int, std::string> _info_message;
    };
}





#endif // __MAST_npsol_optimization_interface_h__
