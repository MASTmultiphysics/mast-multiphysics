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


#ifndef __mast__rigid_surface_motion_h__
#define __mast__rigid_surface_motion_h__


// MAST includes
#include "boundary_condition/surface_motion_base.h"


namespace MAST {
    
    // Forward declerations
    class FrequencyFunction;
    
    
    class RigidSurfaceMotion:
    public MAST::SurfaceMotionBase {
        
    public:
        
        RigidSurfaceMotion(const std::string& nm);
        
        
        virtual ~RigidSurfaceMotion();

        
        /*!
         *   initiate the motion object with oscillation data. 
         *   \par pitch_phase is the phase angle by which pitch leads
         *   plunge.
         */
        void init(MAST::FrequencyFunction& freq,
                  const RealVectorX& plunge_vector,
                  const RealVectorX& pitch_axis,
                  const RealVectorX& hinge_location,
                  const Real plunge_amplitude,
                  const Real pitch_amplitude,
                  const Real pitch_phase);
        
        
        /*!
         *  provides the value of surface velocity, deformation, and
         *  change in surface normal at the specified time and spatial location
         */
        virtual void
        time_domain_motion(const Real t,
                           const libMesh::Point& p,
                           const libMesh::Point& n,
                           RealVectorX& wdot,
                           RealVectorX& dn_rot);

        /*!
         *   provides a function for the definition of surface displacement,
         *   \par w, and rotation of the surface normal in \par dn_rot
         */
        virtual void
        freq_domain_motion(const libMesh::Point& p,
                           const libMesh::Point& n,
                           ComplexVectorX& w,
                           ComplexVectorX& dn_rot);

        
    protected:
        
        
        /*!
         *   frequency of oscillation
         */
        MAST::FrequencyFunction *_freq;
        
        
        /*!
         *   amplitude and phase data for oscillation
         */
        Real
        _plunge_amplitude,
        _pitch_amplitude,
        _pitch_phase;
        

        /*!
         *   axis for definition of rigid motion
         */
        RealVector3
        _plunge_vector,
        _pitch_axis,
        _hinge_location;
        
    };
}


#endif // __mast__rigid_surface_motion_h__
