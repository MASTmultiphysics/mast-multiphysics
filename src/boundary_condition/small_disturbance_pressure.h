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


#ifndef __mast__small_disturbance_pressure_h__
#define __mast__small_disturbance_pressure_h__


// MAST includes
#include "base/field_function_base.h"


// libMesh includes
#include "libmesh/point.h"


namespace MAST {
    
    
    class SmallDisturbancePressure:
    public MAST::FieldFunction<Real> {
        
    public:
        
        SmallDisturbancePressure();
        
        
        virtual ~SmallDisturbancePressure();
        
        
        /*!
         *   provides a function for the definition of surface displacement,
         *   \par w, and rotation of the surface normal in \par dn_rot
         */
        virtual void
        freq_domain_pressure(const libMesh::Point& p,
                             const bool if_cp,
                             Real& press,
                             Complex& dpress) = 0;
        
        
    protected:

        
    };
}


#endif // __mast__small_disturbance_pressure_h__

