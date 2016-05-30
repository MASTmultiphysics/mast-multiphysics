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

#ifndef __mast_surface_motion_base_h__
#define __mast_surface_motion_base_h__

// MAST includes
#include "base/field_function_base.h"

// libMesh includes
#include "libmesh/point.h"


namespace MAST {
    
    
    class SurfaceMotionBase:
    public MAST::FieldFunction<Real> {

        
    public:
        
        SurfaceMotionBase(const std::string& nm);
        
        
        virtual ~SurfaceMotionBase();
        
        /*!
         *   @returns a clone of the function
         */
        virtual std::auto_ptr<MAST::FieldFunction<Real> >
        clone() const {
            
            libmesh_error(); // must be implemented in derived class
        }

        
        /*!
         *  provides the value of surface velocity, deformation, and 
         *  change in surface normal at the specified time and spatial location
         */
        virtual void
        time_domain_motion(const Real t,
                           const libMesh::Point& p,
                           const libMesh::Point& n,
                           RealVectorX& wdot,
                           RealVectorX& dn_rot) = 0;
        
        
        /*!
         *   provides the value of surface displacement,
         *   \par w, and rotation of the surface normal in \par dn_rot
         */
        virtual void
        freq_domain_motion(const libMesh::Point& p,
                           const libMesh::Point& n,
                           ComplexVectorX& w,
                           ComplexVectorX& dn_rot) = 0;
        
    protected:
        
        
    };
}



#endif // __mast_surface_motion_base_h__
