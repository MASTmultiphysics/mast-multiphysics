/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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


#ifndef __mast__plate_flag_frequency_domain_pressure_function_h__
#define __mast__plate_flag_frequency_domain_pressure_function_h__


// MAST includes
#include "fluid/frequency_domain_pressure_function.h"


// libMesh includes
#include "libmesh/system.h"
#include "libmesh/mesh_function.h"


namespace MAST {
    
    
    
    class PlateFlagFrequencyDomainPressureFunction:
    public MAST::FrequencyDomainPressureFunction {
        
    public:
        
        PlateFlagFrequencyDomainPressureFunction(MAST::SystemInitialization& sys,
                                                MAST::FlightCondition&      flt,
                                                const Real th);
        
        
        virtual ~PlateFlagFrequencyDomainPressureFunction();
        
        
        /*!
         *   provides the complex pressure perturbation
         */
        virtual void
        operator() (const libMesh::Point& p,
                    const Real t,
                    Complex& dp) const;
        
        
    protected:
        
        Real _flag_th;
    };
}


#endif // __mast__plate_flag_frequency_domain_pressure_function_h__

