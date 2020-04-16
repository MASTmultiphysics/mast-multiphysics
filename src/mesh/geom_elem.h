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


#ifndef __mast_geom_elem_h__
#define __mast_geom_elem_h__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/fe_type.h"


namespace MAST {
    
    // Forward declerations
    class FEBase;
    class BoundaryConditionBase;
    class SystemInitialization;
    
    /*!
     *   This class acts as a wrapper around libMesh::Elem for the purpose
     *   of providing a uniform interface for cases where:
     *       - two-dimensional elements may exist in this three-dimensional
     *         space, and
     *       - level-set intersection may create sub-elements inside the
     *         element on either side of the interface.
     *    A physics implementation should be averse to these variants.
     *    This class and it children shoudl honor the following definitions:
     *        - a reference element is the element in the finite element mesh. In
     *          a standard FE analysis a quadrature is performed on this element.
     *        - a quadrature element is the element that is a subset of the
     *          reference element on which the quadrature points are defined. In a
     *          standard FE analysis the reference and quadrature elements are the same.
     *        - a local reference element is the element transformed in a coordinate
     *          system that is defined in along the 1D element or in the plane of
     *          the 2D element.
     *        - a local quadrature element is similarly defined for the
     *          quadrature element.
     */
    class GeomElem {
        
    public:
        GeomElem();
        
        virtual ~GeomElem();
        
        /*!
         *   @return a reference to element in the mesh.
         */
        virtual const libMesh::Elem&
        get_reference_elem() const;

        /*!
         *   @return a reference to element in the mesh.
         */
        virtual const libMesh::Elem& get_reference_local_elem() const;

        /*!
         *   @return a reference to quadrature element.
         */
        virtual const libMesh::Elem& get_quadrature_elem() const;

        /*!
         *   @return a reference to quadrature element.
         */
        virtual const libMesh::Elem& get_quadrature_local_elem() const;

        /*!
         *   @return dimension of the element. 
         */
        unsigned int dim() const;

        /*!
         *   number of sides on quadrature element.
         */
        unsigned int n_sides_quadrature_elem() const;
        
        /*!
         *   @return FEType for the underlying FE discretization for the
         *   i-th variable
         */
        libMesh::FEType get_fe_type(unsigned int i) const;

        /*!
         *   for 1D elements the transformed coordinate system attached to the
         *   element defines the local x-axis along the length of the element.
         *   The local-y axis should be specified using this method. This
         *   then provides the complete information about setting up the local
         *   three-dimensional coordinate system.
         */
        void set_local_y_vector(const RealVectorX& y_vec);
        
        /*!
         * This sets the 1D elements to extension/torsional stiffness only.
         * This is useful when modeling truss structures (e.g. using CROD
         * elements in Nastran) which do not require an orientation vector
         * like beam elements do. By default, 1D elements include bending.
         */
        void set_bending(bool onoff);
        
        /*!
         *   initialize the object for the specified reference \p elem.
         */
        virtual void init(const libMesh::Elem& elem,
                          const MAST::SystemInitialization& sys_init);

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
                     bool init_second_order_derivative,
                     int extra_quadrature_order = 0) const;

        
        /*!
         *   From the given list of boundary loads, this identifies the sides
         *   of the quadrature element and the loads in \p bc that are to be
         *   applied on it. The map of side and loads is returned in \p loads.
         */
        
        
        virtual void
        external_side_loads_for_quadrature_elem
        (std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc,
         std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>>& loads) const;
        
        /*!
         *
         */
        virtual void
        get_boundary_ids_on_quadrature_elem_side
        (unsigned int s, std::vector<libMesh::boundary_id_type>& bc_ids) const;

        /*!
         *   Vector and matrix quantities defined on one- and two-dimensional
         *   elements that are oriented in two or three-dimensional spaces may
         *   need to be transformed to/from element coordinate system. This
         *   is required for one-dimensional elements that are not along
         *   the x-axis or two-dimensional elements that are not in the x-y
         *   plane. @return \p true if transformation is required, \p false
         *   otherwise.
         */
        bool use_local_elem() const;


        const RealVectorX&
        domain_surface_normal() const;
        
        void
        transform_point_to_global_coordinate(const libMesh::Point& local_pt,
                                             libMesh::Point& global_pt) const;

        void
        transform_vector_to_global_coordinate(const libMesh::Point& local_vec,
                                              libMesh::Point& global_vec) const;

        void
        transform_vector_to_local_coordinate(const libMesh::Point& global_vec,
                                             libMesh::Point& local_vec) const;

        void
        transform_vector_to_global_coordinate(const RealVectorX& local_vec,
                                              RealVectorX& global_vec) const;
        
        void
        transform_vector_to_local_coordinate(const RealVectorX& global_vec,
                                             RealVectorX& local_vec) const;

        /*!O
         *   @return the 3x3 transformation matrix to transform a vector from
         *   the element to global coordinate system:
         *   \f$ v_{global} = T v_{elem} \f$.
         */
        const RealMatrixX& T_matrix() const;
        
    protected:
        
        
        /*!
         *   initializes the local element
         */
        void _init_local_elem();

        /*!
         *   initializes the local element
         */
        void _init_local_elem_1d();

        /*!
         *   initializes the local element
         */
        void _init_local_elem_2d();

        /*!
         *  system initialization object for this element
         */
        const MAST::SystemInitialization  *_sys_init;

        bool                               _use_local_elem;
        
        /*!
         *   reference element in the mesh for which the data structure is
         *   initialized
         */
        const libMesh::Elem*               _ref_elem;
        
        /*!
         *    a local element is created if
         */
        libMesh::Elem*                     _local_elem;

        /*!
         *   the y-axis that is used to define the coordinate system for
         *   a 1-D element. This must be provided a local element is required.
         */
        RealVectorX                        _local_y;

        /*!
         *   surface normal of the 1D/2D element.
         */
        RealVectorX                        _domain_surface_normal;

        /*!
         *   nodes for local element
         */
        std::vector<libMesh::Node*>        _local_nodes;

        /*!
         *    Transformation matrix defines T_ij = V_i^t . Vn_j, where
         *    V_i are the unit vectors of the global cs, and Vn_j are the
         *    unit vectors of the local cs. To transform a vector from global to
         *    local cs,    an_j = T^t a_i, and the reverse transformation is
         *    obtained as  a_j  = T  an_i
         */
        RealMatrixX                        _T_mat;
        
        /*!
         * Defines if bending is used in this element or not. True by default
         * Added for github issue #40
         */
        bool                                _bending = true;
    };
}

#endif  // __mast_geom_elem_h__
