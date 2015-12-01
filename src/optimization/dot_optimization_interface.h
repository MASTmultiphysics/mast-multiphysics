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

#ifndef __MAST_dot_optimization_interface_h__
#define __MAST_dot_optimization_interface_h__

// MAST includes
#include "optimization/optimization_interface.h"

extern "C" {
    extern void dot_(int*    INFO,
                     int*    METHOD,
                     int*    IPRINT,
                     int*    NDV,
                     int*    NCON,
                     double* X,
                     double* XL,
                     double* XU,
                     double* OBJ,
                     int*    MINMAX,
                     double* G,
                     double* RPRM,
                     int*    IPRM,
                     double* WK,
                     int*    NRWK,
                     int*    IWK,
                     int*    NRIWK);
}


namespace MAST {
    
    class DOTOptimizationInterface: public MAST::OptimizationInterface {
        
    public:
        
        DOTOptimizationInterface()
        { }
        
        virtual ~DOTOptimizationInterface()
        { }
        
        virtual void optimize();
        
        
    protected:
        
    };
}





#endif // __MAST_dot_optimization_interface_h__
