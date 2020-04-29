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

#ifndef __mast__level_set_boundary_velocity_h__
#define __mast__level_set_boundary_velocity_h__

// MAST includes
#include "base/mesh_field_function.h"


namespace MAST {
    
    class LevelSetBoundaryVelocity:
    public MAST::FieldFunction<RealVectorX> {
    public:
        
        LevelSetBoundaryVelocity(const unsigned int dim);
        
        virtual ~LevelSetBoundaryVelocity();
        
        void init(MAST::SystemInitialization& sys,
                  const MAST::MeshFieldFunction& phi);

        void clear();

        virtual void derivative (const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 RealVectorX& v) const;

        void velocity(const MAST::FunctionBase& f,
                      const libMesh::Point& p,
                      const Real t,
                      RealVectorX& v) const;

        /*!
         * attaches the level set function \p phi with this object. This is necessary only when the
         * interface point functions are used.
         */
        void attach_level_set_function(const MAST::FieldFunction<Real>& phi);

        /*!
         * clears the attached level set function
         */
        void clear_level_set_function();
        
        void search_nearest_interface_point(const libMesh::Elem& e,
                                            const unsigned int side,
                                            const libMesh::Point& p,
                                            const Real t,
                                            RealVectorX& pt) const;

        void search_nearest_interface_point_derivative(const MAST::FunctionBase& f,
                                                       const libMesh::Elem& e,
                                                       const unsigned int side,
                                                       const libMesh::Point& p,
                                                       const Real t,
                                                       RealVectorX& v) const;


        /*!
         * serches for a point \p pt in the vicinity of \p p on the level set interface, where
         * level set function is zero. \p length is a reference length that is used to identify
         * the step-size for the search. If the interface point is expected to be within a few elements,
         * then this length coudl be the element edge length.
         */
        void search_nearest_interface_point_old(const libMesh::Point& p,
                                                const Real t,
                                                const Real length,
                                                RealVectorX& pt,
                                                bool allow_sub_search = true) const;

        /*!
         * serches for a point \p pt in the vicinity of \p p on the level set interface, where
         * level set function is zero. \p length is a reference length that is used to identify
         * the step-size for the search. If the interface point is expected to be within a few elements,
         * then this length coudl be the element edge length.
         */
        void search_nearest_interface_point_derivative_old(const MAST::FunctionBase& f,
                                                           const libMesh::Point& p,
                                                           const Real t,
                                                           const Real length,
                                                           RealVectorX& v) const;

        void normal_at_point(const libMesh::Point& p,
                             const Real t,
                             RealVectorX& n) const;
        

        void normal_derivative_at_point(const MAST::FunctionBase& f,
                                        const libMesh::Point& p,
                                        const Real t,
                                        RealVectorX& n) const;
        
        
    protected:
        
        Real _evaluate_point_search_obj(const libMesh::Point& p,
                                        const Real t,
                                        const RealVectorX& dv) const;

        unsigned int                     _dim;
        const MAST::MeshFieldFunction*   _phi;
        libMesh::MeshBase*               _mesh;
        const MAST::FieldFunction<Real>* _level_set_func;
    };
}

#endif // __mast__level_set_boundary_velocity_h__

