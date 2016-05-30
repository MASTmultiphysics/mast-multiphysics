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


#ifndef __mast__function_surface_motion_h__
#define __mast__function_surface_motion_h__


// MAST includes
#include "boundary_condition/surface_motion_base.h"



namespace MAST {
    
    // Forward declerations
    class FrequencyFunction;
    class DisplacementFunctionBase;
    
    /*!
     *   This class provides the surface motion at a given point and time
     *   through a user defined displacement function and frequency.
     */
    class FunctionSurfaceMotion:
    public MAST::SurfaceMotionBase {
        
    public:
        
        FunctionSurfaceMotion(const std::string&            nm);
        
        
        virtual ~FunctionSurfaceMotion();
        
        
        /*!
         *   initiate the motion object with oscillation data.
         */
        void init(MAST::FrequencyFunction&           freq,
                  MAST::DisplacementFunctionBase&    w);
        
        
        /*!
         *  provides the value of surface velocity, deformation, and
         *  change in surface normal at the specified time and spatial location
         */
        virtual void
        time_domain_motion(const Real               t,
                           const libMesh::Point&    p,
                           const libMesh::Point&    n,
                           RealVectorX&             wdot,
                           RealVectorX&             dn_rot);
        
        /*!
         *   provides a function for the definition of surface displacement,
         *   \par w, and rotation of the surface normal in \par dn_rot
         */
        virtual void
        freq_domain_motion(const libMesh::Point&    p,
                           const libMesh::Point&    n,
                           ComplexVectorX&          w,
                           ComplexVectorX&          dn_rot);
        
        
    protected:
        
        
        /*!
         *   frequency of oscillation
         */
        MAST::FrequencyFunction              *_freq;
        
        
        /*!
         *   function that defines the displacement at a given time
         */
        MAST::DisplacementFunctionBase       *_w;
        
    };
}


#endif // __mast__function_surface_motion_h__
