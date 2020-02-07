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

#ifndef __mast__naca0012_h__
#define __mast__naca0012_h__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/unstructured_mesh.h"

namespace MAST {

    namespace Examples {
 
        class NACA0012Mesh2D {
           
        public:

            NACA0012Mesh2D();
            
            virtual ~NACA0012Mesh2D();
            
            void mesh(const Real r,
                      const Real L,
                      const unsigned int radial_divs,
                      const unsigned int quarter_divs,
                      const Real elem_size_ratio,
                      libMesh::UnstructuredMesh& mesh,
                      libMesh::ElemType etype);
            
        protected:
            
            
        };
        
    }
}


#endif // __mast__cylinder_h__
