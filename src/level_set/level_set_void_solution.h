/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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

#ifndef __mast__level_set_void_solution__
#define __mast__level_set_void_solution__


// MAST includes
#include "base/assembly_base.h"


namespace MAST {

    // Forward declerations
    class LevelSetInterfaceDofHandler;
    class LevelSetIntersection;
    
    /*!
     *   This will compute the solution at the interface under the
     *   assumption of zero surface normal flux. The element solution
     *   in the positive region of the level set function is computed
     *   from finite element analysis, while the solution on the
     *   nodes in the negative domain are identified through Schur-Factorization.
     */
    
    class LevelSetVoidSolution:
    public MAST::AssemblyBase::SolverMonitor {
        
    public:
        
        LevelSetVoidSolution();
        
        virtual ~LevelSetVoidSolution();
        
        virtual void init(MAST::AssemblyBase& assembly) {
            libmesh_error(); // call the other method instead
        }

        virtual void init(MAST::AssemblyBase& assembly,
                          MAST::LevelSetIntersection& intersection,
                          MAST::LevelSetInterfaceDofHandler  &dof_handler);
        
        virtual void clear();

        MAST::AssemblyBase&
        get_assembly()  {
            return *_assembly;
        }
        
        void update_void_solution(libMesh::NumericVector<Real>& X,
                                  libMesh::NumericVector<Real>& dX);

        
    protected:

        
        MAST::AssemblyBase                 *_assembly;

        MAST::LevelSetIntersection         *_intersection;

        MAST::LevelSetInterfaceDofHandler  *_dof_handler;
    };
}

#endif // __mast__level_set_void_solution__

