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


#ifndef __mast_level_set_intersected_elem_h__
#define __mast_level_set_intersected_elem_h__

// MAST includes
#include "mesh/geom_elem.h"


namespace MAST {
    
    // Forward declerations
    class LevelSetIntersection;
    
    /*!
     *   This class inherits from MAST::GeomElem and provides an interface
     *   to initialize FE objects on sub-elements obtained from intersection
     *   of level-set function with a reference element.
     */
    class LevelSetIntersectedElem: public MAST::GeomElem  {
        
    public:
        LevelSetIntersectedElem();
        
        virtual ~LevelSetIntersectedElem();
        
        /*!
         *   @return a reference to quadrature element.
         */
        virtual const libMesh::Elem& get_quadrature_elem() const;
        
        /*!
         *   @return a reference to quadrature element.
         */
        virtual const libMesh::Elem& get_quadrature_local_elem() const;

        /*!
         *   @return \p true if the element has a zero level set boundary
         *   that passes through it.
         */
        bool if_elem_has_level_set_boundary() const;

        
        /*!
         *   @return \p true if the sub element has a side that is coincident
         *   with the zero level set boundary.
         */
        bool if_subelem_has_side_on_level_set_boundary() const;

        
        /*!
         *   @return the side id of the sub element which is coincident
         *   with the zero level set boundary. 
         */
        int get_subelem_side_on_level_set_boundary() const;

        /*!
         *   This method should not get called for this class. Call the
         *   overloaded method that accepts the intersection object.
         */
        virtual void init(const libMesh::Elem& elem,
                          const MAST::SystemInitialization& sys_init) {
            
            libmesh_error(); // should not get called
        }
        
        /*!
         *   Initializes the object for a sub element obtained from the
         *   level-set intersection of a reference element.
         */
        virtual void init(const libMesh::Elem& elem,
                          const MAST::SystemInitialization& sys_init,
                          MAST::LevelSetIntersection& intersection);

        
        /*!
         *   initializes the finite element shape function and quadrature
         *   object with the order of quadrature rule changed based on the
         *   \p extra_quadrature_order.
         */
        virtual std::unique_ptr<MAST::FEBase>
        init_fe(bool init_grads,
                bool init_second_order_derivative,
                int extra_quadrature_order = 0) const;
        
        
        /*!
         *   initializes the finite element shape function and quadrature
         *   object for the side with the order of quadrature rule
         *   changed based on the \p extra_quadrature_order
         */
        virtual std::unique_ptr<MAST::FEBase>
        init_side_fe(unsigned int s,
                     bool init_grads,
                     int extra_quadrature_order = 0) const;

        
    protected:
        
        const libMesh::Elem                   *_sub_elem;
        libMesh::Elem                         *_local_sub_elem;

        MAST::LevelSetIntersection            *_intersection;

        /*!
         *   nodes for local element
         */
        std::vector<libMesh::Node*>           _local_subelem_nodes;
    };
}

#endif  // __mast_level_set_intersected_elem_h__
