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

#ifndef __mast_blade_stiffened_panel_mesh__
#define __mast_blade_stiffened_panel_mesh__

// C++ includes
#include <map>

// MAST includes
#include "base/mast_data_types.h"
#include "base/constant_field_function.h"
#include "base/physics_discipline_base.h"


// libMesh includes
#include "libmesh/point.h"
#include "libmesh/mesh_base.h"

namespace MAST {
    
    
    /*!
     *    builds the mesh for a stiffened panel
     */
    class StiffenedPanelMesh {
    public:
        StiffenedPanelMesh() { }
        
        ~StiffenedPanelMesh() { }
        
        
        void init (const unsigned int n_stiff,
                   const unsigned int n_x_divs,
                   const unsigned int n_y_divs_between_stiffeners,
                   const Real length,
                   const Real width,
                   libMesh::MeshBase& mesh,
                   libMesh::ElemType t,
                   bool beam_stiffeners,
                   const unsigned int n_z_divs_stiff = 0,
                   const Real h_stiff = 0.);
        
    protected:
        
        enum Component {
            PANEL,
            STIFFENER_X,
            STIFFENER_Y
        };
        
        void _combine_mesh(libMesh::MeshBase& panel,
                           libMesh::MeshBase& stiffener,
                           MAST::StiffenedPanelMesh::Component c,
                           Real stiff_offset,
                           libMesh::subdomain_id_type sid);
        
    };
}


#endif // __mast_blade_stiffened_panel_mesh__

