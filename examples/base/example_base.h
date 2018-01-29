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

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/getpot.h"
#include "libmesh/fe_type.h"

namespace MAST {
    
    // Forward declerations
    class Parameter;
    class FunctionBase;
    class BoundaryConditionBase;
    
    namespace Examples {
        
        class ExampleBase {
            
        public:
            
            ExampleBase();
            
            virtual ~ExampleBase();
            
            virtual void init(GetPot& input);
            
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
             *   register a loading
             */
            void register_loading(MAST::BoundaryConditionBase& l);

            
            void add_load_parameter(MAST::Parameter& p);
            
            
            /*!
             *   this should
             */
            void update_load_parameters(Real scale);
            
        protected:
            
            /*!
             *   registers a parameter for sensitivity
             */
            void register_paramter_for_sensitivity(MAST::Parameter& p);
            
            bool                             _initialized;
            GetPot*                          _input;
            libMesh::FEType                  _fetype;
            
        private:
            std::vector<MAST::Parameter*>                           _params_for_sensitivity;
            std::map<std::string, MAST::Parameter*>                 _parameters;
            std::map<std::string, MAST::FunctionBase*>              _field_functions;
            std::set<MAST::BoundaryConditionBase*>                  _loadings;
            std::map<MAST::Parameter*, const Real>                  _load_parameters;
        };
    }
}


#endif //__mast_example_base_h__

