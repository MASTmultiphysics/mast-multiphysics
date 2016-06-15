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

#ifndef __mast__local_elem_base__
#define __mast__local_elem_base__


// MAST includes
#include "base/mast_data_types.h"
#include "base/field_function_base.h"


// libMesh includes
#include "libmesh/elem.h"



namespace MAST {

    // Forward declerations
    template <typename ValType> class FieldFunction;
    
    /*!
     *   class provides a simple mechanism to create a geometric element
     *   in the local coordinate system.
     */
    class LocalElemBase {
    public:
        LocalElemBase(const libMesh::Elem& elem):
        _elem(elem),
        _local_elem(nullptr),
        _T_mat_function(nullptr)
        { }
        
        
        virtual ~LocalElemBase();
        
        /*!
         *   returns a constant reference to the global element.
         */
        const libMesh::Elem& global_elem() const {
            return _elem;
        }
        
        
        /*!
         *   returns a constant reference to the local element.
         */
        const libMesh::Elem& local_elem() const {
            if (!_local_elem) // original element lies in the xy-plane
                return _elem;
            else
                return *_local_elem;
        }
        
        /*!
         *    returns the transformation matrix for this element. This is used
         *    to map the coordinates from local to global coordinate system
         */
        const RealMatrixX& T_matrix() const {
            return _T_mat;
        }
        
        
        /*!
         *    returns the transformation matrix for this element. This is used
         *    to map the coordinates from local to global coordinate system
         */
        const MAST::FieldFunction<RealMatrixX>&
        T_matrix_function();
        
        
        /*!
         *    Local elements are defined for 1D and 2D elements that exist in
         *    3D space. These elements have a local coordinate system associated
         *    with the local coordinate. This method accepts the point defined
         *    the local coordinate system as the input and maps it to the
         *    global coordinate system.
         */
        void global_coordinates_location(const libMesh::Point& local,
                                         libMesh::Point& global) const;

        
        /*!
         *   Calculates the side surface normal of the element at the given
         *   point \par p, where the finite element provided surface normal is
         *   \p n_local. The calculated surface normal in the global
         *   coordinate system is \p n_global.
         *
         *   This method is needed because 1D and 2D elements can live in global
         *   3D space. To deal with this a geometric element with a local
         *   coordinate system is defined and used for initialization of
         *   finite elements.
         */
        void global_coordinates_normal(const libMesh::Point& local,
                                       RealVector3& global) const;


        /*!
         *   Calculates the surface normal of the element at the given point
         *   \par p, where the finite element provided surface normal is
         *   \p n_local. The calculated surface normal in the global
         *   coordinate system is \p n_global. For 1D or 2D elements, the flux
         *   loads can be defined over the entire domain of the element. This,
         *   of course, is an idealization of very thin surfaces that are
         *   exchanging heat with its surroundings.
         *
         *   This method is needed because 1D and 2D elements can live in global
         *   3D space. To deal with this a geometric element with a local
         *   coordinate system is defined and used for initialization of
         *   finite elements.
         */
        virtual void domain_surface_normal_in_global_coordinates(const libMesh::Point& p,
                                                                 RealVector3& n_global) const = 0;

    protected:
        
        /*!
         *   given element in global coordinate system
         */
        const libMesh::Elem& _elem;
        
        /*!
         *   element created in local coordinate system
         */
        libMesh::Elem* _local_elem;
        
        /*!
         *   nodes for local element
         */
        std::vector<libMesh::Node*> _local_nodes;
        
        /*!
         *    Transformation matrix defines T_ij = V_i^t . Vn_j, where
         *    V_i are the unit vectors of the global cs, and Vn_j are the
         *    unit vectors of the local cs. To transform a vector from global to
         *    local cs,    an_j = T^t a_i, and the reverse transformation is
         *    obtained as  a_j  = T  an_i
         */
        RealMatrixX _T_mat;

        
        MAST::FieldFunction<RealMatrixX>* _T_mat_function;
    };


    
    
    /*!
     *   Presently, this function is defined for planar elements. 
     *   Extensions to curved elements is to be added in the near future.
     */
    class TransformMatrixFunction:
    public MAST::FieldFunction<RealMatrixX> {
        
    public:
        
        TransformMatrixFunction(const RealMatrixX& Tmat);
        
        virtual ~TransformMatrixFunction() {}
        
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 RealMatrixX& v) const;
        
        
        virtual void derivative (const MAST::DerivativeType d,
                                 const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 RealMatrixX& v) const;
        
    protected:
        

        const RealMatrixX& _Tmat;
    };
}

#endif // __mast__local_elem_base__

