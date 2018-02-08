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
                  const libMesh::NumericVector<Real>& sol,
                  const libMesh::NumericVector<Real>& dsol);
        
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 RealVectorX& v) const;
        
    protected:
        
        unsigned int             _dim;
        MAST::MeshFieldFunction *_phi;
    };
}

#endif // __mast__level_set_boundary_velocity_h__

