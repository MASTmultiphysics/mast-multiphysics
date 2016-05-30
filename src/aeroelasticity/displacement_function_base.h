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

#ifndef __mast__displacement_function_base_h__
#define __mast__displacement_function_base_h__

// MAST includes
#include "base/field_function_base.h"


namespace MAST {
    
    // Forward declerations
    class Parameter;
    
    class DisplacementFunctionBase:
    public MAST::FieldFunction<RealVectorX> {
        
    public:
        
        DisplacementFunctionBase(const std::string& nm):
        MAST::FieldFunction<RealVectorX>(nm) { }
        
        
        virtual ~DisplacementFunctionBase() { }
        
        
        virtual std::auto_ptr<MAST::FieldFunction<RealVectorX> > clone() const {
            libmesh_assert(false);
        }
        
        /*!
         *    calculates the value of the displacement derivative wrt x.
         */
        virtual void dwdx (const libMesh::Point& p,
                           const Real t,
                           RealVectorX& dwdx) const = 0;

        
        /*!
         *    calculates the value of the displacement derivative wrt y.
         */
        virtual void dwdy (const libMesh::Point& p,
                           const Real t,
                           RealVectorX& dwdx) const = 0;

        
        /*!
         *    calculates the value of the displacement derivative wrt z.
         */
        virtual void dwdz (const libMesh::Point& p,
                           const Real t,
                           RealVectorX& dwdz) const = 0;

        
    protected:
        
        
    };
}

#endif // __mast__frequency_function_h__
