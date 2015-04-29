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

#ifndef __mast__system_initialization__
#define __mast__system_initialization__

// C++ includes
#include <memory>


// libMesh includes
#include "libmesh/system.h"


namespace MAST {
    
    class SystemInitialization {
    public:
        /*!
         *   initialize the variables in the provided system \par sys
         *   of \par order and \par family. Uses \par prefix for
         *   all variables name.
         */
        SystemInitialization (libMesh::System& sys,
                              const std::string& prefix);
        
        /*!
         *   virtual destructor
         */
        virtual ~SystemInitialization();
        
        
        /*!
         *   @returns the number of variables in this system
         */
        unsigned int n_vars() const {
            return _system.n_vars();
        }
        
        /*!
         *   @returns the FEType object for variable \par i,
         *   that defines the finite element family and order.
         */
        const libMesh::FEType&
        fetype(unsigned int i) const {
            return _system.variable_type(i);
        }
        
        
        /*!
         *  @returns a reference to the system for which the variables
         *  are initialized.
         */
        libMesh::System& system() {
            return _system;
        }
        
        /*!
         *  @returns a constant reference to the system for which the variables
         *  are initialized.
         */
        const libMesh::System& system() const {
            return _system;
        }
        
        /*!
         *  @returns a constant reference to the vector of variable IDs.
         */
        const std::vector<unsigned int> vars() const {
            return _vars;
        }
        
        /*!
         *  @returns a constant reference to the prefix used for all
         *  variables.
         */
        const std::string& prefix() const {
            return _prefix;
        }
        
        
//        /*!
//         *   @returns a smart-pointer to the solution function of this system.
//         */
//        std::auto_ptr<MAST::FieldFunction<RealVectorX> >
//        solution_function();
        
    protected:
        
        libMesh::System& _system;
        
        std::vector<unsigned int> _vars;
        
        std::string _prefix;
    };
}


#endif //__mast__system_initialization__
