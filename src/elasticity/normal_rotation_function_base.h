/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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

#ifndef __mast__normal_rotation_function_base__
#define __mast__normal_rotation_function_base__

// MAST includes
#include "base/mesh_field_function.h"


namespace MAST {
    
    /*!
     *    This uses the displacement gradient to calculate the rotation in a
     *    given surface normal.
     */
    template <typename ValType>
    class NormalRotationFunctionBase:
    public MAST::FieldFunction<ValType> {
        
    public:

        NormalRotationFunctionBase(const std::string& nm):
        MAST::FieldFunction<ValType>(nm) {
            
        }
        
        
        virtual ~NormalRotationFunctionBase() { }
        
        
        virtual void operator() (const libMesh::Point& p,
                                 const libMesh::Point& n,
                                 const Real t,
                                 ValType& dn_rot) const = 0;
        
        virtual void perturbation (const libMesh::Point& p,
                                   const libMesh::Point& n,
                                   const Real t,
                                   ValType& dn_rot) const = 0;
        
    protected:
        
    };
}

#endif // __mast__normal_rotation_function_base__

