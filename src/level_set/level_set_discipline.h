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

#ifndef __mast__level_set_discipline__
#define __mast__level_set_discipline__

// MAST includes
#include "base/physics_discipline_base.h"



namespace MAST {
    
    // Forward declerations
    template <typename ValType> class FieldFunction;
    
    
    class LevelSetDiscipline:
    public MAST::PhysicsDisciplineBase {
        
        
    public:
        
        // Constructor
        LevelSetDiscipline(libMesh::EquationSystems& eq_sys);
        
        
        /*!
         *   virtual destructor
         */
        virtual ~LevelSetDiscipline();

        
        /*!
         *   sets the level set function normal velocity field
         */
        void
        set_velocity_function(const MAST::FieldFunction<Real>& vel) {
            
            libmesh_assert(!_vel);
            _vel = &vel;
        }

        /*!
         *   @returns true if the velocity function has been set for this
         *   discipline
         */
        bool
        has_velocity_function() const  {
            
            return _vel;
        }
        
        /*!
         *  @returns a reference to the velocity function for this level set
         */
        const MAST::FieldFunction<Real>&
        get_velocity_function() const {
            
            return *_vel;
        }

        /*!
         *   If \p true, then the level set will be propagated using the
         *   velocity specified by the velocity object in the constructor.
         *   Otherwise, the level set will be reinitialized for
         *   \f$ |\nabla(\phi)| = 1 \f$.
         */
        void set_level_set_propagation_mode(bool f) {
            _if_level_set_propagation = f;
        }
        
        
        bool if_level_set_propagation() const {
            
            return _if_level_set_propagation;
        }
        
    protected:
        
        
        const MAST::FieldFunction<Real>* _vel;

        bool _if_level_set_propagation;
    };
}


#endif // __mast__level_set_discipline__

