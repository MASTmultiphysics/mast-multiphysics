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

#ifndef __mast__indicator_function_constrain_dofs_h__
#define __mast__indicator_function_constrain_dofs_h__

// MAST includes
#include "level_set/level_set_constrain_dofs.h"

// libMesh includes
#include "libmesh/system.h"


namespace MAST {
    
    // Forward declerations
    template <typename ValType> class FieldFunction;
    class LevelSetIntersection;
    class SystemInitialization;
    
    /*!
     *   Constrains the dofs based on level set function and indicator function.
     *   Any dofs that are completely on the negative side are  constrained.
     *   Additionally, any dofs that belong to elements with zero function
     *   value are constrained.
     */
    class IndicatorFunctionConstrainDofs:
    public MAST::LevelSetConstrainDofs {
    public:
        
        
        IndicatorFunctionConstrainDofs(MAST::SystemInitialization&        sys,
                                       MAST::FieldFunction<Real>&         level_set,
                                       MAST::FieldFunction<RealVectorX>&  indicator);
        
        
        virtual ~IndicatorFunctionConstrainDofs();
        
        /*!
         *   provides implementation of the libMesh::System::Constraint::constrain()
         *   virtual method
         */
        virtual void
        constrain ();
        
    protected:
        
        MAST::FieldFunction<RealVectorX>            &_indicator;
    };
}


#endif //__mast__indicator_function_constrain_dofs_h__



