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

#ifndef __mast__normal_rotation_mesh_function__
#define __mast__normal_rotation_mesh_function__

// MAST includes
#include "elasticity/normal_rotation_function_base.h"


namespace MAST {

    // Forward declerations
    class MeshFieldFunction;
    
    
    
    class NormalRotationMeshFunction:
    public MAST::NormalRotationFunctionBase<RealVectorX> {
        
    public:
        
        NormalRotationMeshFunction(const std::string& nm,
                                   MAST::MeshFieldFunction& func);
        
        
        virtual ~NormalRotationMeshFunction() { }
        
        
        virtual void operator() (const libMesh::Point& p,
                                 const libMesh::Point& n,
                                 const Real t,
                                 RealVectorX& dn_rot) const;
        
        virtual void perturbation (const libMesh::Point& p,
                                   const libMesh::Point& n,
                                   const Real t,
                                   RealVectorX& dn_rot) const;
        
    protected:

        /*!
         *   mesh field function
         */
        MAST::MeshFieldFunction& _func;
        
    };
}

#endif // __mast__normal_rotation_mesh_function__

