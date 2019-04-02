/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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

#ifndef __mast__beam_flag_frequency_domain_normal_rotation__
#define __mast__beam_flag_frequency_domain_normal_rotation__

// MAST includes
#include "elasticity/complex_normal_rotation_mesh_function.h"


namespace MAST {
    
    
    class BeamFlagFrequencyDomainNormalRotation:
    public MAST::ComplexNormalRotationMeshFunction {
        
    public:
        
        BeamFlagFrequencyDomainNormalRotation(const std::string& nm,
                                              MAST::ComplexMeshFieldFunction& func,
                                              const std::vector<Real>& midplane);
        
        
        virtual ~BeamFlagFrequencyDomainNormalRotation() { }
        
        
        virtual void operator() (const libMesh::Point& p,
                                 const libMesh::Point& n,
                                 const Real t,
                                 ComplexVectorX& dn_rot) const;
        
        virtual void perturbation (const libMesh::Point& p,
                                   const libMesh::Point& n,
                                   const Real t,
                                   ComplexVectorX& dn_rot) const;
        
    protected:
        
        /*!
         *   @returns the y coordinate of the nearest midplane for this point
         */
        Real _nearest_midplane_y_coord(const libMesh::Point& p) const;
        
        /*!
         *   vector of midplane coordinates
         */
        std::vector<Real>   _mid_coords;
    };
}

#endif // __mast__beam_flag_frequency_domain_normal_rotation__


