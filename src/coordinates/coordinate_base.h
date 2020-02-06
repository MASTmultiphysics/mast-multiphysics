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


#ifndef __mast__coordinate_base__
#define __mast__coordinate_base__


// MAST includes
#include "base/field_function_base.h"


namespace MAST {
    
    /*!
     *    Provides the transformation matrix T to transform
     *    vector from the orientation provided in this matrix,
     *    to one in the global basis
     */
    class CoordinateBase:
    public MAST::FieldFunction<RealMatrixX> {
        
    public:
        
        CoordinateBase(const std::string& nm);

        /*!
         *   prepares the matrix \p mat that transforms stress and strain tensors
         *   represented in a 6x1 vector from the coordinate system in _orient
         *   to the global coordinate system. Note that the shear straints in
         *   the strain tensor vector should be represented in the tensor quantities,
         *   and not the engineering strain.
         */
        void stress_strain_transformation_matrix(const RealMatrixX& T,
                                                 RealMatrixX& mat) const;
        

        void stress_strain_transformation_matrix_sens(const RealMatrixX& T,
                                                      const RealMatrixX& dT,
                                                      RealMatrixX& mat) const;
    protected:
        
    };
}




#endif // __mast__coordinate_base__
