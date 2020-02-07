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

#ifndef __mast__flag_mesh_2D_h__
#define __mast__flag_mesh_2D_h__

// MAST includes
#include "examples/fluid/meshing/mesh_initializer.h"

// libMesh includes
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"


namespace MAST {
    
    
    class FlagMesh3D: public MAST::MeshInitializer
    {
    public:
        FlagMesh3D():
        _flag_th    (0.),
        _x_in       (0.),         _x_out(0.),
        _x_le       (0.),          _x_te(0.),
        _y_flag_left(0.),  _y_flag_right(0.),
        _y_dom_left (0.),   _y_dom_right(0.),
        _z_lo(0.),  _z_up(0.),
        _panel_bc_id (0)
        { }
        
        
        Real thickness() const {return _flag_th;}
        
        /*!
         *   initializes the object with the division for each dimension.
         */
        virtual void
        init (const unsigned int panel_bc_id,
              const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
              libMesh::UnstructuredMesh& mesh, libMesh::ElemType t);
        
    protected:
        
        virtual void process_mesh( );
        
        /*!
         *   dimensions used to create mesh
         */
        Real _flag_th;
        
        Real _x_in, _x_out, _x_le, _x_te;

        Real _y_flag_left, _y_flag_right, _y_dom_left, _y_dom_right;
        
        Real _z_lo, _z_up;
        
        unsigned int _panel_bc_id;
    };
    
}


#endif //__mast__flag_mesh_2D_h__
