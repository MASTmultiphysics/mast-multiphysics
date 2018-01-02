/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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

#ifndef __mast__plate_flag_frequency_domain_displacement__
#define __mast__plate_flag_frequency_domain_displacement__

// MAST includes
#include "base/complex_mesh_field_function.h"



namespace MAST {
    
    
    /*!
     *    This provides a wrapper FieldFunction compatible class that
     *    interpolates the solution using libMesh's MeshFunction class.
     */
    class PlateFlagFrequencyDomainDisplacement:
    public MAST::ComplexMeshFieldFunction {
        
    public:
        /*!
         *   constructor
         */
        PlateFlagFrequencyDomainDisplacement(MAST::SystemInitialization& sys,
                                            const std::string& nm);
        
        
        /*!
         *   destructor
         */
        virtual ~PlateFlagFrequencyDomainDisplacement();
        
        
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 ComplexVectorX& v) const;
        
        
        virtual void perturbation (const libMesh::Point& p,
                                   const Real t,
                                   ComplexVectorX& v) const;
        
        
    protected:
        
    };
}

#endif // __mast__plate_flag_frequency_domain_displacement__


