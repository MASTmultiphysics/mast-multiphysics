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

#ifndef __mast_cylinder_analysis_2d_h__
#define __mast_cylinder_analysis_2d_h__


// MAST includes
#include "examples/fluid/base/fluid_example_base.h"

namespace MAST {
    
    namespace Examples {
        
        
        class CylinderAnalysis2D:
        public MAST::Examples::FluidExampleBase {
            
        public:
            
            CylinderAnalysis2D(const libMesh::Parallel::Communicator& comm_in);
            
            virtual ~CylinderAnalysis2D();
            
        protected:
            
            virtual void _init_mesh();
            virtual void _init_loads();
        };
    }
}



#endif //  __mast_cylinder_analysis_2d_h__
