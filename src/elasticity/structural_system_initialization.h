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

#ifndef __mast__structural_system_initialization__
#define __mast__structural_system_initialization__

// MAST includes
#include "base/system_initialization.h"

// libMesh includes
#include "libmesh/system.h"

namespace MAST {
    
    // Forward declerations
    class NonlinearSystem;
    
    class StructuralSystemInitialization:
    public MAST::SystemInitialization  {
        
    public:
        StructuralSystemInitialization(MAST::NonlinearSystem& sys,
                                       const std::string& prefix,
                                       const libMesh::FEType& fe_type);
        
        virtual ~StructuralSystemInitialization();
        
        /*!
         *   @returns a reference to libMesh::System that stores the stress
         *   variables
         */
        libMesh::System& get_stress_sys() {
            
            return *_stress_output_sys;
        }
        
        /*!
         *   @returns a reference to vector of strain-stress variable ids in the
         *   sequence \f$ \{
         *   \epsilon_{xx}, \epsilon_{yy}, \epsilon_{zz},
         *   \epsilon_{xy}, \epsilon_{yz}, \epsilon_{zx},
         *   \sigma_{xx},     \sigma_{yy},   \sigma_{zz},
         *   \sigma_{xy},     \sigma_{yz},   \sigma_{zx},
         *   \sigma_{vm}
         *   \} \f$
         */
        const std::vector<unsigned int>& get_stress_var_ids() {
            
            return _stress_vars;
        }

        
    protected:
        
        libMesh::System*           _stress_output_sys;
        
        std::vector<unsigned int>  _stress_vars;
    };
}

#endif // __mast__structural_system_initialization__
