/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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

// MAST includes
#include "coordinates/coordinate_base.h"


MAST::CoordinateBase::CoordinateBase(const std::string& nm):
MAST::FieldFunction<RealMatrixX>(nm) {
    
}



void
MAST::CoordinateBase::
stress_strain_transformation_matrix(const RealMatrixX& T,
                                    RealMatrixX &mat) const {
    
    mat.setZero(6,6);
    
    std::vector<std::pair<unsigned int, unsigned int> > ids(6);
    ids[0] = std::pair<unsigned int, unsigned int>(0,0);
    ids[1] = std::pair<unsigned int, unsigned int>(1,1);
    ids[2] = std::pair<unsigned int, unsigned int>(2,2);
    ids[3] = std::pair<unsigned int, unsigned int>(0,1);
    ids[4] = std::pair<unsigned int, unsigned int>(1,2);
    ids[5] = std::pair<unsigned int, unsigned int>(2,0);
    
    for (unsigned int i=0; i<6; i++)
        for (unsigned int j=0; j<3; j++) {
            // the first three columns of the matrix
            mat(i,j)   =
            T(ids[j].first, ids[i].first) * T(ids[j].second, ids[i].second);
            //            libMesh::out
            //            << "( " << i << ", " << j << " ): "
            //            << ids[j].first+1 << ids[i].first+1 << " x "
            //            << ids[j].second+1 << ids[i].second+1 << std::endl;
            
            // last three columns of the matrix
            mat(i,j+3) =
            T(ids[j+3].first, ids[i].first)  * T(ids[j+3].second, ids[i].second) +
            T(ids[j+3].first, ids[i].second) * T(ids[j+3].second, ids[i].first);
            
            //            libMesh::out
            //            << "( " << i << ", " << j+3 << " ): "
            //            << ids[j+3].first+1 << ids[i].first+1  << " x "
            //            <<  ids[j+3].second+1 << ids[i].second+1
            //            << " + " << ids[j+3].first+1 << ids[i].second+1 << " x "
            //            << ids[j+3].second+1 << ids[i].first+1 << std::endl;
        }
}



void
MAST::CoordinateBase::
stress_strain_transformation_matrix_sens(const RealMatrixX& T,
                                         const RealMatrixX& dT,
                                         RealMatrixX &mat) const {
    
    mat.setZero(6,6);
    
    std::vector<std::pair<unsigned int, unsigned int> > ids(6);
    ids[0] = std::pair<unsigned int, unsigned int>(0,0);
    ids[1] = std::pair<unsigned int, unsigned int>(1,1);
    ids[2] = std::pair<unsigned int, unsigned int>(2,2);
    ids[3] = std::pair<unsigned int, unsigned int>(0,1);
    ids[4] = std::pair<unsigned int, unsigned int>(1,2);
    ids[5] = std::pair<unsigned int, unsigned int>(2,0);
    
    for (unsigned int i=0; i<6; i++)
        for (unsigned int j=0; j<3; j++) {
            // the first three columns of the matrix
            mat(i,j)   =
            dT(ids[j].first, ids[i].first) *  T(ids[j].second, ids[i].second) +
            T (ids[j].first, ids[i].first) * dT(ids[j].second, ids[i].second);
            
            // last three columns of the matrix
            mat(i,j+3) =
            dT(ids[j+3].first, ids[i].first)  *  T(ids[j+3].second, ids[i].second) +
            T (ids[j+3].first, ids[i].first)  * dT(ids[j+3].second, ids[i].second) +
            dT(ids[j+3].first, ids[i].second) *  T(ids[j+3].second, ids[i].first) +
            T (ids[j+3].first, ids[i].second) * dT(ids[j+3].second, ids[i].first);
        }
}


