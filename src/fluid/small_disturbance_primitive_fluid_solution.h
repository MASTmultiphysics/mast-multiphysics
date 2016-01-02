/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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

#ifndef __mast__small_disturbance_primitive_fluid_solution__
#define __mast__small_disturbance_primitive_fluid_solution__


// MAST includes
#include "base/mast_data_types.h"


namespace MAST {

    // Forward declerations
    class PrimitiveSolution;
    
    /*!
     *   Class defines basic operations and calculation of the small 
     *   disturbance primitive variables.
     */
    template <typename ValType>
    class SmallPerturbationPrimitiveSolution
    {
    public:
        SmallPerturbationPrimitiveSolution();
        
        void zero();
        
        
        void init(const MAST::PrimitiveSolution& sol,
                  const typename VectorType<ValType>::return_type& delta_sol);
        
        
        void print(std::ostream& out) const;
        
        
        ValType c_pressure(const Real q0) const;
        
        
        void get_duvec(typename VectorType<ValType>::return_type& du) const;
        
        
        typename VectorType<ValType>::return_type perturb_primitive_sol;
        
        
        ValType drho;
        
        ValType du1;
        
        ValType du2;
        
        ValType du3;
        
        ValType dT;
        
        ValType dp;
        
        ValType da;
        
        ValType de_tot;
        
        ValType dk;
        
        ValType dentropy;
        
        ValType dmach;
        
        const MAST::PrimitiveSolution* primitive_sol;
    };
    
    
}


#endif /* defined(__mast__small_disturbance_primitive_fluid_solution__) */
