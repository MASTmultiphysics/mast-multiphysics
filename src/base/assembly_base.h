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

#ifndef __mast__system_assembly_base__
#define __mast__system_assembly_base__


// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/system.h"


namespace MAST {
    
    // Forward declerations
    class PhysicsDisciplineBase;
    class SystemInitialization;
    class ElementBase;
    
    
    class AssemblyBase {
    public:
        
        /*!
         *  constructor takes a reference to the discipline that provides
         *  the boundary conditions, volume loads, properties, etc.
         */
        AssemblyBase();
        
        /*!
         *   virtual destructor
         */
        virtual ~AssemblyBase();
        
        /*!
         *   @returns a const reference to the PhysicsDisciplineBase object
         *   associated with this object
         */
        const MAST::PhysicsDisciplineBase& discipline() const;
        
        /*!
         *   @returns a non-const reference to the PhysicsDisciplineBase object
         *   associated with this object
         */
        MAST::PhysicsDisciplineBase& discipline();
        
        
        /*!
         *   attaches a system to this discipline, and vice-a-versa
         */
        virtual void
        attach_discipline_and_system(MAST::PhysicsDisciplineBase& discipline,
                                     MAST::SystemInitialization& system) = 0;

        
        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void
        clear_discipline_and_system( ) = 0;

        
        /*!
         *   @returns a const reference to the libMesh::System object
         *   associated with this object
         */
        const libMesh::System& system() const;
        
        /*!
         *   @returns a non-const reference to the libMesh::System object
         *   associated with this object
         */
        libMesh::System& system();
        
        
    protected:
        
        /*!
         *   @returns a smart-pointer to a newly created element for
         *   calculation of element quantities.
         */
        virtual std::auto_ptr<MAST::ElementBase>
        _build_elem(const libMesh::Elem& elem) = 0;
        
        /*!
         *   localizes the parallel vector so that the local copy
         *   stores all values necessary for calculation of the
         *   element quantities
         */
        std::auto_ptr<libMesh::NumericVector<Real> >
        _build_localized_vector(const libMesh::System& sys,
                                const libMesh::NumericVector<Real>& global);
        
        
        /*!
         *   PhysicsDisciplineBase object for which this class is assembling
         */
        MAST::PhysicsDisciplineBase* _discipline;
        
        
        /*!
         *   System for which this assembly is performed
         */
        MAST::SystemInitialization* _system;
    };
        
}

#endif // __mast__system_assembly_base__

