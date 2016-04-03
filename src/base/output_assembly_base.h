/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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

#ifndef __mast__output_assembly_base__
#define __mast__output_assembly_base__

// MAST includes
#include "base/assembly_base.h"


namespace MAST {
    
    
    class OutputAssemblyBase:
    public MAST::AssemblyBase,
    public libMesh::System::QOI,
    public libMesh::System::QOIDerivative,
    public libMesh::System::QOIParameterSensitivity {
    
    public:
        OutputAssemblyBase();
        
        virtual ~OutputAssemblyBase();
        
        /*!
         *   attaches a system to this discipline, and vice-a-versa
         */
        virtual void
        attach_discipline_and_system(MAST::PhysicsDisciplineBase& discipline,
                                     MAST::SystemInitialization& system) ;
        
        
        
        /*!
         *   Reattaches to the same system that was attached earlier.
         *
         *   This cannot be called if the clear_discipline_and_system() method
         *   has been called.
         */
        virtual void
        reattach_to_system();
        

        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void
        clear_discipline_and_system( );
        
        
        
    protected:
        
    };
}


#endif // __mast__output_assembly_base__

