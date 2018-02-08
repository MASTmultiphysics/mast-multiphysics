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

#ifndef __mast_example_base_h__
#define __mast_example_base_h__

// C++ includes
#include <string>
#include <vector>
#include <memory>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/getpot.h"
#include "libmesh/fe_type.h"
#include "libmesh/parallel_object.h"


namespace MAST {
    
    // Forward declerations
    class Parameter;
    class FunctionBase;
    class BoundaryConditionBase;
    
    namespace Examples {
        
        // Forward declerations
        class GetPotWrapper;
        
        class ExampleBase:
        public libMesh::ParallelObject {
            
        public:
            
            ExampleBase(const libMesh::Parallel::Communicator& comm_in);
            
            virtual ~ExampleBase();
            
            virtual void init(MAST::Examples::GetPotWrapper& input,
                              const std::string& prefix);

            /*!
             *  If the user says, this is prepended to the input parameters
             *  names.
             */
            std::string prefix();
            
            /*!
             *   adds a parameter
             */
            void add_parameter(MAST::Parameter& p);
            
            /*!
             *   @returns a parameter by the specified name
             */
            MAST::Parameter& get_parameter(const std::string &nm);
            
            /*!
             *   adds a function
             */
            void register_field_function(MAST::FunctionBase& f);

            /*!
             *   register a boundary condition
             */
            void register_loading(MAST::BoundaryConditionBase& l);

            
            void add_load_parameter(MAST::Parameter& p);
            
            
            /*!
             *   this should update the load parameter to a value between 0 and 1.
             *
             */
            void update_load_parameters(Real scale);
            
        protected:
            
            /*!
             *   registers a parameter for sensitivity
             */
            void register_paramter_for_sensitivity(MAST::Parameter& p);
            
            bool                             _initialized;
            std::string                      _prefix;
            MAST::Examples::GetPotWrapper*   _input;
            libMesh::FEType                  _fetype;
            
        private:
            std::vector<MAST::Parameter*>                           _params_for_sensitivity;
            std::map<std::string, MAST::Parameter*>                 _parameters;
            std::vector<MAST::Parameter*>                           _dv_parameters;
            std::map<std::string, MAST::FunctionBase*>              _field_functions;
            std::set<MAST::BoundaryConditionBase*>                  _boundary_conditions;
            std::map<MAST::Parameter*, const Real>                  _load_parameters;
        };
    }
}


#endif //__mast_example_base_h__

