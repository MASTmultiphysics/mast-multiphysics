/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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

#ifndef __mast__panel_mesh_2D_h__
#define __mast__panel_mesh_2D_h__

// MAST includes
#include "examples/fluid/meshing/mesh_initializer.h"

// libMesh includes
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"


namespace MAST {
    
    
    class PanelMesh2D: public MAST::MeshInitializer
    {
    public:
        PanelMesh2D():
        _tc_ratio(0.),
        _x0(0.), _x1(0.),
        _y0(0.), _y1(0.),
        _n_maxima(0),
        _panel_bc_id(0),
        _symmetry_bc_id(0),
        _cos_profile(false)
        { }
        
        /*!
         *   initializes the object with the division for each dimension. Sets-up the
         *   mesh for the provided information, and then uses the provided funciton
         *   move the mesh points. If cos_profile = false, the function used to
         *   define the panel surface is sin (n_maxima * pi * x / L ). Otherwise,
         *   the function is  1 - cos(n_maxima * 2 * pi * x / L)
         */
        virtual void init (const Real tc, bool cos_profile,
                           const unsigned int n_maxima,
                           const unsigned int panel_bc_id,
                           const unsigned int symmetry_bc_id,
                           const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                           libMesh::UnstructuredMesh& mesh, libMesh::ElemType t);
        
    protected:
        
        /*!
         *   move the mesh points to account for a bump, and apply boudnary condition
         */
        virtual void process_mesh( );
        
        /*!
         *   t/c ratio of the panel
         */
        Real _tc_ratio;
        
        Real _x0, _x1;
        
        Real _y0, _y1;
        
        unsigned int _n_maxima;
        
        unsigned int _panel_bc_id, _symmetry_bc_id;
        
        bool _cos_profile;
    };
    
}


#endif // __mast__panel_mesh_2D_h__
