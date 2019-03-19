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

#ifndef __mast_assembly_elem_operation_h__
#define __mast_assembly_elem_operation_h__

// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/elem.h"

namespace MAST {
    
    // Forward declerations
    class FEBase;
    class ElementBase;
    class AssemblyBase;
    class FunctionBase;
    class SystemInitialization;
    class PhysicsDisciplineBase;
    
    class AssemblyElemOperations {
        
    public:
        AssemblyElemOperations();
        
        virtual ~AssemblyElemOperations();
        
        
        /*!
         *   @returns a reference to the system initialization object
         */
        MAST::SystemInitialization& get_system_initialization();
        
        /*!
         *   @returns a reference to the discipline object
         */
        MAST::PhysicsDisciplineBase& get_discipline();
        
        /*!
         *   attaches a system to this discipline
         */
        virtual void
        set_discipline_and_system(MAST::PhysicsDisciplineBase& discipline,
                                  MAST::SystemInitialization& system);
        
        
        /*!
         *   clears association with a system to this discipline
         */
        virtual void
        clear_discipline_and_system();

        /*!
         *   sets the assembly object
         */
        virtual void set_assembly(MAST::AssemblyBase& assembly);
        
        /*!
         *   @returns a reference to the assembly object
         */
        virtual MAST::AssemblyBase& get_assembly();
        
        /*!
         *   clears the assembly object
         */
        virtual void clear_assembly();
        
        /*!
         *   initializes the object for calculation of element quantities for
         *   the specified \p elem.
         */
        virtual void
        init(const libMesh::Elem& elem) = 0;
        
        /*!
         *   clears the element initialization
         */
        virtual void clear_elem();
        
        /*!
         *   @returns a reference to the physics element. The object must have
         *   been initialized before this method is called.
         */
        MAST::ElementBase& get_physics_elem() {
            
            libmesh_assert(_physics_elem);
            return *_physics_elem;
        }
         
        /*!
         *   sets the element solution
         */
        virtual void
        set_elem_solution(const RealVectorX& sol);

        /*!
         *   sets the element solution sensitivity
         */
        virtual void
        set_elem_solution_sensitivity(const RealVectorX& sol);

        /*!
         *   sets the element perturbed solution
         */
        virtual void
        set_elem_perturbed_solution(const RealVectorX& sol);

        
        /*!
         *   sets the element velocity
         */
        virtual void set_elem_velocity(const RealVectorX& vel);

        /*!
         *   sets the element velocity sensitivity
         */
        virtual void set_elem_velocity_sensitivity(const RealVectorX& vel);

        /*!
         *   sets the element perturbed velocity
         */
        virtual void set_elem_perturbed_velocity(const RealVectorX& vel);

        /*!
         *   sets the element acceleration
         */
        virtual void set_elem_acceleration(const RealVectorX& accel);
        
        /*!
         *   sets the element acceleration
         */
        virtual void set_elem_acceleration_sensitivity(const RealVectorX& accel);

        /*!
         *   sets the element perturbed acceleration
         */
        virtual void set_elem_perturbed_acceleration(const RealVectorX& accel);

        
    protected:

        MAST::SystemInitialization       *_system;
        MAST::PhysicsDisciplineBase      *_discipline;

        MAST::AssemblyBase               *_assembly;
        
        MAST::ElementBase                *_physics_elem;
    };
}


#endif // __mast_assembly_elem_operation_h__

