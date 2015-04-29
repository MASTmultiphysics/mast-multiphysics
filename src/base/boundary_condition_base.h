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

#ifndef __mast__boundary_condition__
#define __mast__boundary_condition__

// MAST includes
#include "base/function_set_base.h"


namespace MAST {
    
    enum BoundaryConditionType {
        
        SURFACE_PRESSURE,
        SMALL_DISTURBANCE_DISPLACEMENT, // displacement perturbations
        SMALL_DISTURBANCE_PRESSURE,     // pressure perturbations about steady values
        SMALL_DISTURBANCE_MOTION,       // pressure and motion perturbations about steady state values
        PISTON_THEORY,
        DIRICHLET,
        TEMPERATURE,
        HEAT_FLUX,
        CONVECTION_HEAT_FLUX,
        SURFACE_RADIATION_HEAT_FLUX,
        HEAT_SOURCE,
        NO_SLIP_WALL,
        SYMMETRY_WALL,
        SLIP_WALL,
        FAR_FIELD,
        EXHAUST,
        ISOTHERMAL,
        ADIABATIC
    };
    
    
    class BoundaryConditionBase:
    public MAST::FunctionSetBase {
        
    public:
        BoundaryConditionBase(MAST::BoundaryConditionType t):
        MAST::FunctionSetBase(),
        _bc_type(t)
        { }
        
        virtual ~BoundaryConditionBase() { }
        
        
        MAST::BoundaryConditionType type() const {
            return _bc_type;
        }
        
    protected:
        
        MAST::BoundaryConditionType _bc_type;
    };
    
}


#endif  // __mast__boundary_condition__
