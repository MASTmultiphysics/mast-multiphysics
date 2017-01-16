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


#ifndef __mast__polar_coordinate_base__
#define __mast__polar_coordinate_base__


// MAST includes
#include "coordinates/coordinate_base.h"


namespace MAST {
    
    /*!
     *    Defines a polar coordinate system with the radius and 
     *    obtained from the two parameters provided in the constructor
     */
    class PolarCoordinate:
    public MAST::CoordinateBase {
        
    public:
        PolarCoordinate(const std::string& nm,
                        MAST::FieldFunction<Real>& theta);
        
        virtual ~PolarCoordinate();
        
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 RealMatrixX& v) const;
        
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \par p, and time, \par t, and returns it in \p v.
         */
        virtual void derivative (const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 RealMatrixX& v) const;
        
        
    protected:
        
        MAST::FieldFunction<Real>& _theta;
    };
}




#endif // __mast__polar_coordinate_base__
