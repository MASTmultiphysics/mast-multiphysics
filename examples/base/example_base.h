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
#include "libmesh/unstructured_mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"


namespace MAST {
    
    // Forward declerations
    class Parameter;
    class FunctionBase;
    class BoundaryConditionBase;
    class NonlinearSystem;
    class SystemInitialization;
    class PhysicsDisciplineBase;
    class MaterialPropertyCardBase;
    class ElementPropertyCardBase;
    class DirichletBoundaryCondition;
    class OutputAssemblyElemOperations;

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

            virtual void _init_mesh() = 0;
            virtual void _init_system_and_discipline() = 0;
            virtual void _init_dirichlet_conditions() = 0;
            virtual void _init_boundary_dirichlet_constraint(const unsigned int bid,
                                                             const std::string& tag) = 0;
            virtual void _init_eq_sys() = 0;
            virtual void _init_loads() = 0;
            virtual void _init_material() = 0;
            virtual void _init_section_property() = 0;

            /*!
             *   registers a parameter for sensitivity
             */
            void register_paramter_for_sensitivity(MAST::Parameter& p);
            
            bool                                                    _initialized;
            std::string                                             _prefix;
            MAST::Examples::GetPotWrapper*                          _input;
            libMesh::FEType                                         _fetype;
            
            libMesh::UnstructuredMesh*                              _mesh;
            libMesh::EquationSystems*                               _eq_sys;
            MAST::NonlinearSystem*                                  _sys;
            MAST::SystemInitialization*                             _sys_init;
            MAST::PhysicsDisciplineBase*                            _discipline;
            MAST::MaterialPropertyCardBase*                         _m_card;
            MAST::ElementPropertyCardBase*                          _p_card;

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

