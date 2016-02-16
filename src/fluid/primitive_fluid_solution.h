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

#ifndef __mast__primitive_fluid_solution_h__
#define __mast__primitive_fluid_solution_h__

// MAST include
#include "base/mast_data_types.h"


namespace MAST {
    

    /*!
     *  Class defines the conversion and some basic operations on 
     *  primitive fluid variables used in calculation of flux, Jacobians, etc.
     */
    class PrimitiveSolution {
        
    public:
        PrimitiveSolution();
        
        void zero();
        
        void init(const unsigned int dim,
                  const RealVectorX& conservative_sol,
                  const Real cp_val,
                  const Real cv_val,
                  bool if_viscous);
        
        void print(std::ostream& out) const;
        
        Real c_pressure(const Real p0, const Real q0) const;
        
        void get_uvec(RealVectorX& u) const;
        
        RealVectorX primitive_sol;
        
        unsigned int dimension;
        
        Real cp;
        
        Real cv;
        
        Real rho;
        
        Real u1;
        
        Real u2;
        
        Real u3;
        
        Real T;
        
        Real p;
        
        Real a;
        
        Real e_tot;
        
        Real k;
        
        Real entropy;
        
        Real mach;
        
        // viscous quantities
        Real Pr;
        
        Real k_thermal;
        
        Real mu;
        
        Real lambda;
    };
    
}

#endif // __mast__primitive_fluid_solution_h__
