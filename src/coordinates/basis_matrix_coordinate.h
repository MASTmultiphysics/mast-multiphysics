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


#ifndef __mast__basis_matrix_coordinate__
#define __mast__basis_matrix_coordinate__


// MAST includes
#include "coordinates/coordinate_base.h"


namespace MAST {
    
    /*!
     *    Provides the transformation matrix T to transform
     *    vector from the orientation provided in this matrix,
     *    to one in the global basis
     */
    class BasisMatrixCoordinate:
    public MAST::CoordinateBase {
        
    public:
        BasisMatrixCoordinate(const std::string& nm,
                              MAST::FieldFunction<RealMatrixX>& basis);
        
        virtual ~BasisMatrixCoordinate();
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \p p, and time, \p t, and returns it in \p v.
         */
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 RealMatrixX& v) const;
        
        
        /*!
         *    calculates the value of the function at the specified point,
         *    \p p, and time, \p t, and returns it in \p v.
         */
        virtual void derivative (const MAST::FunctionBase& f,
                                 const libMesh::Point& p,
                                 const Real t,
                                 RealMatrixX& v) const;
        
        
    protected:
        
        MAST::FieldFunction<RealMatrixX>& _basis;
        
    };
}




#endif // __mast__basis_matrix_coordinate__

