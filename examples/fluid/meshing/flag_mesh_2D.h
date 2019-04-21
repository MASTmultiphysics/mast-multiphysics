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

#ifndef __mast__flag_mesh_2D_h__
#define __mast__flag_mesh_2D_h__

// C++ includes
#include <vector>

// MAST includes
#include "examples/fluid/meshing/mesh_initializer.h"

// libMesh includes
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"


namespace MAST {
    
    
    class FlagMesh2D: public MAST::MeshInitializer
    {
    public:
        FlagMesh2D(const unsigned int n_flags,
                   const unsigned int panel_bc_id,
                   const std::vector<MeshInitializer::CoordinateDivisions*>& divs);


        /*!
         *   @returns a reference to the vector of mid-plane coordinate of
         *   all flags
         */
        const std::vector<Real>& midplane_coordinates() const {return _flag_center;}

        
        /*!
         *   @returns the mid-plane coordinate of i^th flag
         */
        Real midplane(unsigned int i) const {return _flag_center[i];}


        /*!
         *   @returns the thickness of i^th flag
         */
        Real thickness(unsigned int i) const {return _flag_th[i];}
        
        /*!
         *   initializes the object with the division for each dimension.
         */
        virtual void
        init_fluid_mesh (libMesh::UnstructuredMesh& mesh,
                         libMesh::ElemType t);

        
        /*!
         *   initializes the structural mesh consistent with the fluid mesh
         */
        virtual void
        init_structural_mesh (libMesh::UnstructuredMesh& mesh,
                              libMesh::ElemType t);

    protected:
        
        virtual void process_mesh( );
        
        libMesh::Elem*
        copy_elem(libMesh::Elem& e,
                  std::map<libMesh::Node*, libMesh::Node*>& added_nodes);
        
        /*!
         *   thickness of flag
         */
        std::vector<Real> _flag_th;

        /*!
         *   centerline of flag
         */
        std::vector<Real> _flag_center;

        Real _x_in, _x_out, _x_le, _x_te;
        
        Real _y_lo, _y_up;

        unsigned int _n_flags;
        
        unsigned int _panel_bc_id;

        const std::vector<MeshInitializer::CoordinateDivisions*>& _divs;
    };
    
}


#endif //__mast__flag_mesh_2D_h__
