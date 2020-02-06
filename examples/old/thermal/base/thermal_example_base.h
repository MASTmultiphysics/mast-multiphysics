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


#ifndef __mast__thermal_example_base_h__
#define __mast__thermal_example_base_h__

// MAST includes
#include "examples/base/example_base.h"



namespace MAST {
    
    
    namespace Examples {
        
        class ThermalExampleBase:
        public MAST::Examples::ExampleBase {
            
        public:
            
            ThermalExampleBase(const libMesh::Parallel::Communicator& comm_in);
            
            virtual ~ThermalExampleBase();
            
            /*!
             *   This is called before each nonlinear and transient solve.
             *   By default, this zeros the solution, but the user can
             *   override it to provide a different initialization.
             */
            virtual void initialize_solution();
            
            
            virtual void steady_solve();
            virtual void steady_sensitivity_solve(MAST::Parameter& p);
            
            virtual void transient_solve();
            virtual void transient_sensitivity_solve(MAST::Parameter& p);
            
        protected:
            
            virtual void _init_system_and_discipline();
            virtual void _init_dirichlet_conditions();
            virtual void _init_boundary_dirichlet_constraint(const unsigned int bid,
                                                             const std::string& tag);
            virtual void _init_eq_sys();
            virtual void _init_loads();
            virtual void _init_material();
            virtual void _init_section_property();
            virtual void _init_flux_load       (bool on_side, unsigned int id_num);
            virtual void _init_convection_load (bool on_side, unsigned int id_num);
            virtual void _init_radiation_load  (bool on_side, unsigned int id_num);
            virtual void _init_source_load     (unsigned int domain_num);
        };
    }
}

#endif // __mast__thermal_example_base_h__

