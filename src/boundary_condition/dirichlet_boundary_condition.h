/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
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

#ifndef __mast__dirichlet_boundary_condition__
#define __mast__dirichlet_boundary_condition__

// C++ includes
#include <memory>

// MAST includes
#include "base/boundary_condition_base.h"


// libmesh includes
#include "libmesh/dirichlet_boundaries.h"



namespace MAST {
    
    // Forward declerations
    template <typename ValType> class FieldFunction;
    
    
    class DirichletBoundaryCondition:
    public MAST::BoundaryConditionBase {
      
    public:
        DirichletBoundaryCondition():
        MAST::BoundaryConditionBase(MAST::DIRICHLET)
        { }
        
        virtual ~DirichletBoundaryCondition() { }

        /*!
         *   initializes the object for the specified domain id (either boundary,
         *   or subdomain), for the displacement components initialized using
         *   a bitwise operator. This method initializes the components to zero
         *   value on the domain.
         */
        void init(const libMesh::boundary_id_type bid,
                  const std::vector<unsigned int>& constrained_vars,
                  MAST::FieldFunction<Real>* f = NULL);
        
        
        /*!
         *    @returns reference to the Dirichlet boundary condition object
         */
        libMesh::DirichletBoundary& dirichlet_boundary() {
            return *_dirichlet_boundary;
        }
        
    protected:
        
        /*!
         *    Dirichlet boundary function for this boundary
         */
        std::auto_ptr<libMesh::DirichletBoundary> _dirichlet_boundary;
    };
}

#endif  // __mast__dirichlet_boundary_condition__
