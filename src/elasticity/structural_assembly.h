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

#ifndef __mast__structural_assembly__
#define __mast__structural_assembly__

// MAST includes
#include "base/assembly_base.h"



namespace MAST {

    // Forward declerations
    class StructuralElementBase;

    /*!
     *   This class provides some routines that are common to
     *   structural assembly routines.
     */
    class StructuralAssembly:
    MAST::AssemblyBase::SolverMonitor {
        
    public:
        
        StructuralAssembly();
        
        virtual ~StructuralAssembly();
        
        virtual void init(MAST::AssemblyBase& assembly);
        
        virtual void clear();

        void set_elem_incompatible_sol(MAST::StructuralElementBase& elem);
        
        void update_incompatible_solution(libMesh::NumericVector<Real>& X,
                                          libMesh::NumericVector<Real>& dX);
        
    protected:
        
        
//        /*!
//         *   assembles the point loas
//         */
//        void _assemble_point_loads(MAST::PhysicsDisciplineBase& discipline,
//                                   MAST::SystemInitialization& system,
//                                   libMesh::NumericVector<Real>& res);

        MAST::AssemblyBase* _assembly;
        
        /*!
         *   map of local incompatible mode solution per 3D elements
         */
        std::map<const libMesh::Elem*, RealVectorX> _incompatible_sol;

    };
}


#endif // __mast__structural_assembly__
