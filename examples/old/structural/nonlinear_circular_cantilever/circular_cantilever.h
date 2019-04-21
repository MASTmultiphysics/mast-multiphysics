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

#ifndef __mast__structural_example_3d_h__
#define __mast__structural_example_3d_h__

// MAST includes
#include "examples/structural/base/structural_example_base.h"


namespace MAST {
    
    namespace Examples {
        
        
        class StructuralExample3D:
        public MAST::Examples::StructuralExampleBase {
            
        public:
            
            StructuralExample3D(const libMesh::Parallel::Communicator& comm_in);
            
            virtual ~StructuralExample3D();
            
        protected:
            
            virtual void _init_loads();
            virtual void _init_mesh();
            virtual void _init_dirichlet_conditions();
            virtual void _init_section_property();
        };
    }
}


#endif  // __mast__structural_example_3d_h__
