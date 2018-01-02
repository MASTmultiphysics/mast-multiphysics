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


#ifndef __mast_local_elem_fe_h__
#define __mast_local_elem_fe_h__


// MAST includes
#include "mesh/fe_base.h"


namespace MAST {

    // Forward declerations
    class LocalElemBase;
    
    /*!
     *   This specialize the MAST::FEBase class for problems where a 1D/2D
     *   element exists in a 3D space, which typically happens in thermal
     *   and structural problems. The class will create a geomtric element
     *   in the 1D/2D space of the geometric element and use that to
     *   initialize the finite element object.
     */
    class LocalElemFE:
    public MAST::FEBase {
      
    public:
        
        LocalElemFE(const MAST::SystemInitialization& sys);
        
        virtual ~LocalElemFE();
        
        /*!
         *   @returns a non-constant reference to the element in the local
         *   coordinate. This is needed for 1D or 2D elements that live
         *   in a 3D space.
         */
        MAST::LocalElemBase& local_elem() {
            return *_local_elem;
        }
        
        
        /*!
         *   @returns a constant reference to the element in the local
         *   coordinate. This is needed for 1D or 2D elements that live
         *   in a 3D space.
         */
        const MAST::LocalElemBase& local_elem() const {
            return *_local_elem;
        }
        
        

        /*!
         *   Creates a local element and then initializes the FE object.
         */
        virtual void init(const libMesh::Elem& elem,
                          const std::vector<libMesh::Point>* pts = nullptr);
        
        /*!
         *   Creates a local element and then initializes the FE object for
         *   the side.
         */
        virtual void init_for_side(const libMesh::Elem& elem,
                                   unsigned int s,
                                   bool if_calculate_dphi);
        
        /*!
         *    a 1D vector only provides the orientation of the x-axis.
         *    It requires a second vector to identify the orientation
         *    of the y-axis. This must be set for 1D elements before
         *    initialization.
         */
        void set_1d_y_vector(const libMesh::Point& y);


        /*!
         *    @returns a constant reference to y-vector for this element.
         *    This must have been set earlier.
         */
        const libMesh::Point& get_1d_y_vector() const;

        /*!
         *    @returns a vector to the global xyz location of the quadrature
         *    points. This is overloaded from MAST::FEBase since this method
         *    transforms the location from local to global coordinate system
         *    for use by the elements.
         */
        virtual const std::vector<libMesh::Point>&
        get_xyz() const {
            
            libmesh_assert(_initialized);
            return _global_xyz;
        }

        virtual const std::vector<libMesh::Point>&
        get_normals() const {
            
            libmesh_assert(_initialized);
            return _global_normals;
        }


    protected:

        void _create_local_element(const libMesh::Elem& elem);
        
        MAST::LocalElemBase   *_local_elem;
        
        libMesh::Point         _y_vector_1d_elem;
        
        std::vector<libMesh::Point> _global_xyz;
        
        std::vector<libMesh::Point> _global_normals;
    };
}


#endif // __mast_local_elem_fe_base_h__
