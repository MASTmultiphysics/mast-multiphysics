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

#ifndef __mast__piston_theory_boundary_condition__
#define __mast__piston_theory_boundary_condition__

// MAST includes
#include "base/boundary_condition_base.h"


namespace MAST {
    
    
    class PistonTheoryBoundaryCondition:
    public MAST::BoundaryConditionBase {
        
    public:
        
        /*!
         *  Constructor for the Piston Theory boundary condition
         *  object. The arguments needed for initialization are
         *  \p order: order of piston theory, \p mach: mach number
         *  \p a_inf: ambient speed of sound, \p gamma: ratio of 
         *  specific heats at constant pressure and constant volume,
         *  \p rho: ambient density for calculation of dynamic pressure,
         *  \p vel_vec: velocity unit vector
         */
        PistonTheoryBoundaryCondition(unsigned int order,
                                      Real mach,
                                      Real a_inf,
                                      Real gamma,
                                      Real rho,
                                      const RealVectorX& vel_vec);
        
        
        virtual ~PistonTheoryBoundaryCondition();

        
        /*!
         *  @returns the order of piston theory to be used
         */
        unsigned int order() const;
        

        /*!
         *  @returns value of ambient Mach number
         */
        Real mach() const;
        

        /*!
         *  @returns value of ambient speed of sound
         */
        Real a_inf() const;
        
        /*!
         *  @returns the ambient density
         */
        Real rho() const;
        
        /*!
         *  @returns value of ambient ratio of specif heat values at
         *  constant pressure and constant volume
         */
        Real gamma() const;


        /*!
         *  @returns velocity vector
         */
        const RealVectorX& vel_vec() const;

    protected:
        

        /*!
         *   Order of the boundary condition
         */
        const unsigned int _order;
        
        
        /*!
         *   Ambient flow property: Mach number
         */
        const Real _mach;


        /*!
         *   Ambient flow property: speed of sound
         */
        const Real _a_inf;
        
        
        /*!
         *   Ambient flow property: ratio of specific heat values at
         *   constant pressure and constant volume
         */
        const Real _gamma;
        
        /*!
         *   Ambient flow property: density
         */
        const Real _rho;

        /*!
         *   Ambient flow velocity vector
         */
        RealVectorX _vel_vec;

    };
}


#endif /* __mast__piston_theory_boundary_condition__ */

