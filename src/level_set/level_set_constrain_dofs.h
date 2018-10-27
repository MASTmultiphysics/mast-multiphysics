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

#ifndef __mast__level_set_constrain_dofs_h__
#define __mast__level_set_constrain_dofs_h__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/system.h"


namespace MAST {
    
    // Forward declerations
    template <typename ValType> class FieldFunction;
    class LevelSetIntersection;
    class SystemInitialization;
    
    /*!
     *   constrains the dofs based on level set function. Any dofs that are
     *   completely on the negative side are constrained.
     */
    class LevelSetConstrainDofs:
    public libMesh::System::Constraint {
    public:
        
        
        LevelSetConstrainDofs(MAST::SystemInitialization& sys,
                              MAST::FieldFunction<Real>& level_set);
        
        
        virtual ~LevelSetConstrainDofs();
        
        
        /*!
         *  if set to true, then if the element has a node on negative side
         * of the level set, we will constrain all dofs of the node to be zero.
         *  If false,
         *  then nodes will be constrained only if they get no contribution
         *  from elements with a positive level set region.
         */
        void constrain_all_negative_indices(bool f) {
            _constrain_all_negative_indices = f;
        }
        
        /*!
         *  @returns a reference to the level set function
         */
        MAST::LevelSetIntersection& get_intersection();
        
        
        /*!
         *   provides implementation of the libMesh::System::Constraint::constrain()
         *   virtual method
         */
        virtual void
        constrain ();
        
    protected:
        
        bool                                  _constrain_all_negative_indices;
        
        MAST::SystemInitialization           &_sys;
        
        MAST::FieldFunction<Real>            &_level_set;
        
        MAST::LevelSetIntersection           *_intersection;
    };
}


#endif //__mast__level_set_constrain_dofs_h__


