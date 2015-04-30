/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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

#ifndef __mast__structural_assembly__
#define __mast__structural_assembly__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/petsc_nonlinear_solver.h"



namespace MAST {

    // Forward declerations
    class PhysicsDisciplineBase;
    class SystemInitialization;

    // monitor function for PETSc solver so that
    // the incompatible solution can be updated after each converged iterate
    PetscErrorCode
    _snes_structural_nonlinear_assembly_monitor_function(SNES snes,
                                                         PetscInt its,
                                                         PetscReal norm2,
                                                         void* ctx);
    

    /*!
     *   This class provides some routines that are common to
     *   structural assembly routines.
     */
    class StructuralAssembly {
        
    public:
        
        StructuralAssembly();
        
        virtual ~StructuralAssembly();
        
        
    protected:
        
        
        /*!
         *   assembles the point loas
         */
        void _assemble_point_loads(MAST::PhysicsDisciplineBase& discipline,
                                   MAST::SystemInitialization& system,
                                   libMesh::NumericVector<Real>& res);
        
    };
}


#endif // __mast__structural_assembly__
