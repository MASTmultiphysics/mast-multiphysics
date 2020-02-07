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

// MAST includes
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "base/field_function_base.h"


// libMesh includes
#include "libmesh/zero_function.h"

void
MAST::DirichletBoundaryCondition::init(const libMesh::boundary_id_type bid,
                                       const std::vector<unsigned int>& constrained_vars) {
    
    // should not have been initialized if this is called
    libmesh_assert(_dirichlet_boundary.get() == nullptr);
    
    std::set<libMesh::boundary_id_type> bid_set;
    bid_set.insert(bid);
    
    std::unique_ptr<libMesh::FunctionBase<Real> > function;

    // if the function was not give, then assume it to be zero function
    function.reset(new libMesh::ZeroFunction<Real>);
    
    _dirichlet_boundary.reset(new libMesh::DirichletBoundary(bid_set,
                                                             constrained_vars,
                                                             function.get()));
}


