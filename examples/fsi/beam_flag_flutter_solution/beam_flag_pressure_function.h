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


#ifndef __mast__beam_flag_pressure_function_h__
#define __mast__beam_flag_pressure_function_h__


// MAST includes
#include "fluid/pressure_function.h"


// libMesh includes
#include "libmesh/system.h"
#include "libmesh/mesh_function.h"


namespace MAST {
    
    // Forward declerations
    class FrequencyFunction;
    class SystemInitialization;
    class FlightCondition;
    
    
    class BeamFlagPressureFunction:
    public MAST::PressureFunction {
        
    public:
        
        BeamFlagPressureFunction(MAST::SystemInitialization& sys,
                                 MAST::FlightCondition&      flt,
                                 const Real th);
        
        
        virtual ~BeamFlagPressureFunction();
        
        
        /*!
         *   provides the value of the pressure at the specified point and time
         */
        virtual void
        operator() (const libMesh::Point& p,
                    const Real t,
                    Real& press) const;
        
        
        /*!
         *   provides the pressure perturbation. The user must have initialized
         *   the perturbed solution using the appropriate init routine.
         */
        virtual void
        perturbation(const libMesh::Point& p,
                     const Real t,
                     Real& dpress) const;
        
        
    protected:
        
        /*!
         *    the flag thickness is assumed to be uniformly distributed about
         *    y=0
         */
        Real _flag_thickness;
    };
}


#endif // __mast__beam_flag_pressure_function_h__

