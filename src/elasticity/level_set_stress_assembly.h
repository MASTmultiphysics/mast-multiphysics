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

#ifndef __mast__level_set_stress_assembly__
#define __mast__level_set_stress_assembly__

// MAST includes
#include "elasticity/stress_assembly.h"


namespace MAST {
    
    // Forward declerations
    template <typename ValType> class FieldFunction;
    class LevelSetIntersection;
    class LevelSetInterfaceDofHandler;

    
    class LevelSetStressAssembly:
    public MAST::StressAssembly {
    public:
        
        
        /*!
         *   constructor associates this assembly object with the system
         */
        LevelSetStressAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~LevelSetStressAssembly();
        
        
        /*!
         *   attaches level set function to \p this. If \p dof_handler is
         *   provided then its solution will be used in computing stress.
         */
        virtual void
        init(MAST::FieldFunction<Real>& level_set,
             MAST::LevelSetInterfaceDofHandler* dof_handler = nullptr);

        
        /*!
         *   clears association with level set function
         */
        virtual void
        clear();

        
        /*!
         *   updates the stresses and strains for the specified solution
         *   vector \p X. Only the maximum values out of each element are
         *   updated. This will put the stress data in the System::solution
         *   vector related to stress/strain values.
         */
        virtual void
        update_stress_strain_data(MAST::StressStrainOutputBase&       ops,
                                  const libMesh::NumericVector<Real>& X);
        

    protected:

        MAST::FieldFunction<Real>            *_level_set;
        MAST::LevelSetIntersection           *_intersection;
        MAST::LevelSetInterfaceDofHandler    *_dof_handler;
    };
}


#endif //__mast__level_set_stress_assembly__


