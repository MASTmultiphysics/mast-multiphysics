/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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

#ifndef __mast__structural_system_h__
#define __mast__structural_system_h__

// C++ includes
#include <memory>

// MAST includes
#include "base/nonlinear_system.h"

// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/enum_eigen_solver_type.h"
#include "libmesh/eigen_system.h"


namespace MAST {
    
    
    // Forward declerations
    class Parameter;
    
    /*!
     *   This class implements a system for quasi-static analysis of 
     *   nonlinear structures.
     */
    class StructuralSystem:
    public MAST::NonlinearSystem {

    public:
        
        /*!
         *    Default constructor
         */
        StructuralSystem(libMesh::EquationSystems& es,
                         const std::string& name,
                         const unsigned int number);
        
        
        virtual ~StructuralSystem();
        
        
        /*!
         *   sets the laod parameter for incrementing
         */
        void set_load_parameter(MAST::Parameter& param,
                                Real min_p,
                                Real max_p);
        
        
        /*!
         * Clear all the data structures associated with
         * the system.
         */
        virtual void clear () libmesh_override;
        
        
        /**
         * Assembles & solves the nonlinear system R(x) = 0.
         */
        virtual void solve () libmesh_override;
        
        
    protected:
        
        
        /*!
         *   iteration counter
         */
        unsigned int _iter;
        

        /*!
         *    value of beta that scales the external load vector in the constrain
         */
        Real _beta;
        
        
        /*!
         *    value of update length
         */
        Real _dl;
        
        
        /*!
         *    load parameter that is updated by this solution procedure
         */
        MAST::Parameter *_load_param;
      
        
        /*!
         *   maximum and minimum values of the parameter
         */
        Real  _min_p, _max_p;
        
    };
}


#endif // __mast__structural_system_h__
