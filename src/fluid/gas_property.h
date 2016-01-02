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

#ifndef mast_gas_property_h
#define mast_gas_property_h


// MAST includes
#include "base/mast_data_types.h"


namespace MAST {
    
    
    class GasProperty
    {
    public:
        GasProperty():
        rho(0.),
        T(0.),
        pressure(0.),
        gamma(0.),
        cp(0.),
        cv(0.),
        R(0.),
        a(0.),
        Pr(0.),
        k_thermal(0.),
        mu(0.),
        lambda(0.)
        {}
        
        
        void zero()
        {
            rho       = 0.;
            T         = 0.;
            pressure  = 0.;
            gamma     = 0.;
            cp        = 0.;
            cv        = 0.;
            R         = 0.;
            a         = 0.;
            Pr        = 0.;
            k_thermal = 0.;
            mu        = 0.;
            lambda    = 0.;
        }
        
        
        /*!
         *   Property values for ideal gas
         */
        Real rho, T, pressure, gamma, cp, cv, R, a;
        
        /*!
         *   Properties for viscous analysis
         */
        Real Pr, k_thermal, mu, lambda;
        
        /*!
         *   initializes the data
         */
        void init();
        
    };
    
    
    
    inline void
    GasProperty::init()
    {
        // the following data should have been set
        libmesh_assert(rho > 0.);
        libmesh_assert(T  >  0.);
        libmesh_assert(cp  > 0.);
        libmesh_assert(cv  > 0.);
        
        R        = cp-cv;
        gamma    = cp/cv;
        pressure = rho*R*T;
        a        = sqrt(gamma*R*T);
    }
    
}

#endif
