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

#ifndef __mast__conservative_fluid_system_initialization_h__
#define __mast__conservative_fluid_system_initialization_h__

// MAST includes
#include "base/system_initialization.h"


namespace MAST {
    
    class ConservativeFluidSystemInitialization:
    public MAST::SystemInitialization  {
        
    public:
        ConservativeFluidSystemInitialization(MAST::NonlinearSystem& sys,
                                              const std::string& prefix,
                                              const libMesh::FEType& fe_type,
                                              const unsigned int dim);
        
        virtual ~ConservativeFluidSystemInitialization();
        
        
        /*!
         *    @returns spatial dimensions of analysis
         */
        inline unsigned int dim() {
            
            return _dim;
        }
        
        
    protected:
        
        
        /*!
         *   spatial dimensions of analysis
         */
        const unsigned int _dim;
    };
}

#endif  //__mast__conservative_fluid_system_initialization_h__


