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

#ifndef __MAST_optimization_interface_h__
#define __MAST_optimization_interface_h__

// MAST includes
#include "base/mast_data_types.h"


namespace MAST {

    // Forward declerations
    class FunctionEvaluation;
    
    /*!
     *    Provides the basic interface API for classes the provide 
     *    implement optimization problems.
     */
    class OptimizationInterface {
    public:
     
        OptimizationInterface():
        _feval(nullptr)
        { }
        
        virtual ~OptimizationInterface()
        { }

        
        virtual void optimize() = 0;
        

        virtual void
        attach_function_evaluation_object (MAST::FunctionEvaluation& feval);
        

        virtual void
        set_real_parameter    (const std::string& nm, Real val) {}

        virtual void
        set_integer_parameter (const std::string& nm, int val) {}

    protected:
        
        MAST::FunctionEvaluation* _feval;
    };
}



#endif  // __MAST_optimization_interface_h__
