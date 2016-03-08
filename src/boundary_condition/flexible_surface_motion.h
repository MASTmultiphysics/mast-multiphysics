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


#ifndef __mast__flexible_surface_motion_h__
#define __mast__flexible_surface_motion_h__


// MAST includes
#include "boundary_condition/surface_motion_base.h"


// libMesh includes
#include "libmesh/system.h"
#include "libmesh/mesh_function.h"


namespace MAST {
    
    // Forward declerations
    class FrequencyFunction;
    class SystemInitialization;
    
    
    class FlexibleSurfaceMotion:
    public MAST::SurfaceMotionBase {
        
    public:
        
        FlexibleSurfaceMotion(MAST::SystemInitialization& sys);
        
        
        virtual ~FlexibleSurfaceMotion();
        
        
        /*!
         *   initiate the motion object with oscillation data.
         *   \par pitch_phase is the phase angle by which pitch leads
         *   plunge.
         */
        void init(MAST::FrequencyFunction& freq,
                  libMesh::NumericVector<Real>& sol);

        
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
         *   system associated with the mesh and solution vector
         */
        MAST::SystemInitialization&         _system;

        
        /*!
         *   frequency of oscillation
         */
        MAST::FrequencyFunction              *_freq;
        
        
        /*!
         *   mesh function that interpolates the solution
         */
        std::auto_ptr<libMesh::MeshFunction> _function;
        
        /*!
         *    numeric vector that stores the solution
         */
        std::auto_ptr<libMesh::NumericVector<Real> > _sol;
    };
}


#endif // __mast__flexible_surface_motion_h__
