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
         *   value on the domain. If \p f_val is not provided then a zero value is
         *   assumed for all constrained variables. \p FieldFunction should
         *   implement a method that provides a vector containing the value of constrained
         *   variables at specified spatial point and time. The order of variables in this vector is
         *   based on the value of \p index. If \p index = \p libMesh::SYSTEM_VARIABLE_ORDER,
         *   which is the default, then the \p i^th component of the vector returned by \p f_val
         *   is the value of the \p i^th system variable. On the other hand, if
         *   If \p index = \p libMesh::LOCAL_VARIABLE_ORDER, then the \p i^th component
         *   of the vector returned by \p f_val is the value of the \p i^th component in
         *   \p constrained_vars. If \p index = \p libMesh::SYSTEM_VARIABLE_ORDER
         *   then \p n_vars should be set to the number of variables in the system.
         */
        void init(const libMesh::boundary_id_type bid,
                  const std::vector<unsigned int>& constrained_vars,
                  MAST::FieldFunction<RealVectorX>* f_val = nullptr,
                  libMesh::VariableIndexing index = libMesh::SYSTEM_VARIABLE_ORDER,
                  unsigned int n_sys_vars = 0);
        
        
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
        std::unique_ptr<libMesh::DirichletBoundary> _dirichlet_boundary;
    };
}

#endif  // __mast__dirichlet_boundary_condition__
