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

#ifndef __mast__structural_near_null_vector_space_h__
#define __mast__structural_near_null_vector_space_h__

// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"


namespace MAST {
    
    /*!
     *  this defines the near-null space of a structural finite element 
     *  model, which is composed of the six rigid-body nodes.
     */
    
    class StructuralNearNullVectorSpace:
    public libMesh::NonlinearImplicitSystem::ComputeVectorSubspace {
        
    public:
        
        /*!
         *   default constructor
         */
        StructuralNearNullVectorSpace();
        

        /*!
         *   destructor
         */
        virtual ~StructuralNearNullVectorSpace() { }
        
        
        /**
         * This function will be called to compute the subspace basis
         * (e.g., nullspace or nearnullspace).
         * It must be implemented by the user in a derived class.
         */
        virtual void operator()(std::vector<libMesh::NumericVector<Real>*>&sp,
                                libMesh::NonlinearImplicitSystem& s);

        
        
    protected:
        
        std::vector<std::string> _nm;
    };
    
}



#endif // __mast__structural_near_null_vector_space_h__

