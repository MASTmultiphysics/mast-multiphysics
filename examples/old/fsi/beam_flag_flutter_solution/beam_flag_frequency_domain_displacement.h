/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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

#ifndef __mast__beam_flag_frequency_domain_displacement__
#define __mast__beam_flag_frequency_domain_displacement__

// C++ includes
#include <vector>

// MAST includes
#include "base/complex_mesh_field_function.h"



namespace MAST {
    
    
    /*!
     *    This provides a wrapper FieldFunction compatible class that
     *    interpolates the solution using libMesh's MeshFunction class.
     */
    class BeamFlagFrequencyDomainDisplacement:
    public MAST::ComplexMeshFieldFunction {
        
    public:
        /*!
         *   constructor
         */
        BeamFlagFrequencyDomainDisplacement(MAST::SystemInitialization& sys,
                                            const std::string& nm,
                                            const std::vector<Real>& midplane);
        
        
        /*!
         *   destructor
         */
        virtual ~BeamFlagFrequencyDomainDisplacement();
        
        
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 ComplexVectorX& v) const;
        
        
        virtual void perturbation (const libMesh::Point& p,
                                   const Real t,
                                   ComplexVectorX& v) const;
        
        
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

#endif // __mast__beam_flag_frequency_domain_displacement__


