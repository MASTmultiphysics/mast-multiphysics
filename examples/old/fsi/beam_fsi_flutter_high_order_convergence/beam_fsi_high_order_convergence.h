/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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

#ifndef __mast_beam_euler_flutter_high_order_convergence_h__
#define __mast_beam_euler_flutter_high_order_convergence_h__

// MAST includes
#include "base/mast_data_types.h"



namespace MAST {

    
    struct BeamFSIFlutterHighOrderConvergence {
    
        
        BeamFSIFlutterHighOrderConvergence();
        
        
        ~BeamFSIFlutterHighOrderConvergence();
        
        
        /*!
         *  solves the system and returns the flutter velocity
         */
        Real solve(bool if_write_output = false,
                   const Real tol = 1.e-1,
                   const unsigned int max_bisection_iters = 20);
        
    };
}



#endif //  __mast_beam_euler_flutter_high_order_convergence_h__

