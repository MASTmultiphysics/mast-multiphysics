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

#ifndef __mast__system_initialization__
#define __mast__system_initialization__

// C++ includes
#include <memory>
#include <vector>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/fe_type.h"


namespace MAST {

    // Forward declerations
    class NonlinearSystem;
    template <typename ValType> class FieldFunction;
    
    
    class SystemInitialization {
    public:
        /*!
         *   initialize the variables in the provided system \par sys
         *   of \par order and \par family. Uses \par prefix for
         *   all variables name.
         */
        SystemInitialization (MAST::NonlinearSystem& sys,
                              const std::string& prefix);
        
        /*!
         *   virtual destructor
         */
        virtual ~SystemInitialization();
        
        
        /*!
         *   @returns the number of variables in this system
         */
        unsigned int n_vars() const;
        
        /*!
         *   @returns the FEType object for variable \par i,
         *   that defines the finite element family and order.
         */
        const libMesh::FEType& fetype(unsigned int i) const;
        
        
        /*!
         *  @returns a reference to the system for which the variables
         *  are initialized.
         */
        MAST::NonlinearSystem& system() {
            return _system;
        }
        
        /*!
         *  @returns a constant reference to the system for which the variables
         *  are initialized.
         */
        const MAST::NonlinearSystem& system() const {
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
        

        /*!
         *    initializes the FE solution vector to the constant
         *    solution provided in \par sol.
         */
        void initialize_solution(const RealVectorX& sol);

        /*!
         *    initializes the FE solution vector to the function
         *    solution provided in \par sol.
         */
        void initialize_solution(const MAST::FieldFunction<RealVectorX>& sol);

        
    protected:
        
        MAST::NonlinearSystem& _system;
        
        std::vector<unsigned int> _vars;
        
        std::string _prefix;
    };
}


#endif //__mast__system_initialization__
