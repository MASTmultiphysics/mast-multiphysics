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

#ifndef __mast__naca0012_wing_h__
#define __mast__naca0012_wing_h__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/unstructured_mesh.h"

namespace MAST {

    namespace Examples {
 
        class NACA0012WingMesh3D {
           
        public:

            NACA0012WingMesh3D();
            
            virtual ~NACA0012WingMesh3D();
            
            void mesh(const Real root_chord,
                      const Real taper_ratio,
                      const Real mid_chord_sweep,
                      const Real far_field_radius_to_root_chord,
                      const Real span,
                      const Real spanwise_farfield,
                      const unsigned int radial_divs_chord,
                      const unsigned int radial_divs_chord_to_farfield,
                      const unsigned int quarter_divs,
                      const unsigned int spanwise_divs,
                      const unsigned int span_to_farfield_divs,
                      const Real radial_elem_size_ratio,
                      const Real spanwise_elem_size_ratio,
                      libMesh::UnstructuredMesh& mesh,
                      libMesh::ElemType etype);
            
        protected:
            
            
        };
        
    }
}


#endif // __mast__naca0012_wing_h__
