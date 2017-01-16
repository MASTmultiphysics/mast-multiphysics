///*
// * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
// * Copyright (C) 2013-2017  Manav Bhatia
// *
// * This library is free software; you can redistribute it and/or
// * modify it under the terms of the GNU Lesser General Public
// * License as published by the Free Software Foundation; either
// * version 2.1 of the License, or (at your option) any later version.
// *
// * This library is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// * Lesser General Public License for more details.
// *
// * You should have received a copy of the GNU Lesser General Public
// * License along with this library; if not, write to the Free Software
// * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// */
//
//#ifndef __mast__shock_impingement_panel_mesh_2D_h__
//#define __mast__shock_impingement_panel_mesh_2D_h__
//
//// MAST includes
//#include "examples/fluid/meshing/mesh_initializer.h"
//
//// libMesh includes
//#include "libmesh/mesh_base.h"
//#include "libmesh/elem.h"
//
//
//namespace MAST {
//    
//    
//    class ShockImpingementPanelMesh2D: public MAST::MeshInitializer
//    {
//    public:
//        ShockImpingementPanelMesh2D():
//        _n_panels(0),
//        _mach(0.),
//        _panel_length(0.),
//        _shock_location_on_panel(0.),
//        _panel_bc_id(0),
//        _symmetry_bc_id(0)
//        { }
//        
//        /*!
//         *   initializes the object with the division for each dimension.
//         */
//        virtual void init (const unsigned int n_panels,
//                           const Real mach,
//                           const Real panel_length,
//                           const Real shock_location_on_panel,
//                           const unsigned int panel_bc_id,
//                           const unsigned int symmetry_bc_id,
//                           const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
//                           libMesh::UnstructuredMesh& mesh, libMesh::ElemType t);
//        
//    protected:
//        
//        virtual void process_mesh( );
//        
//        unsigned int _n_panels;
//
//        Real _mach;
//        
//        Real _panel_length;
//        
//        Real _shock_location_on_panel;
//        
//        unsigned int _panel_bc_id, _symmetry_bc_id;
//    };
//    
//}
//
//
//#endif // __mast__shock_impingement_panel_mesh_2D_h__
