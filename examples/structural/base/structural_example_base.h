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


#ifndef __mast__structural_example_base_h__
#define __mast__structural_example_base_h__

// MAST includes
#include "examples/base/example_base.h"



namespace MAST {
    
    // Forward declerations
    class TimeDomainFlutterSolver;
    class FlutterRootBase;
    
    namespace Examples {
        
        class StructuralExampleBase:
        public MAST::Examples::ExampleBase {
            
        public:
            
            StructuralExampleBase(const libMesh::Parallel::Communicator& comm_in);
            
            virtual ~StructuralExampleBase();
            
            /*!
             *   This is called before each nonlinear and transient solve.
             *   By default, this zeros the solution, but the user can
             *   override it to provide a different initialization.
             */
            virtual void initialize_solution();
            

            virtual void static_solve();
            virtual void static_sensitivity_solve(MAST::Parameter& p);
            virtual void static_adjoint_sensitivity_solve(//MAST::OutputAssemblyElemOperations& q,
                                                          MAST::Parameter& p);
            
            virtual void modal_solve(std::vector<Real>& eig);
            virtual void modal_sensitivity_solve(MAST::Parameter& p, std::vector<Real>& deig_dp);
            virtual void modal_solve_with_nonlinear_load_stepping();
            
            virtual void transient_solve();
            virtual void transient_sensitivity_solve(MAST::Parameter& p);

            virtual void piston_theory_flutter_solve();
            virtual void piston_theory_flutter_sensitivity_solve(MAST::Parameter& p);
            
        protected:
            
            virtual void _init_system_and_discipline();
            virtual void _init_dirichlet_conditions();
            virtual void _init_boundary_dirichlet_constraint(const unsigned int bid,
                                                             const std::string& tag);
            virtual void _init_eq_sys();
            virtual void _init_loads();
            virtual void _init_material();
            virtual void _init_section_property();
            virtual void _init_pressure_load(bool on_side, unsigned int id_num);
            virtual void _init_temperature_load();
            virtual void _init_piston_theory_load();
            
            /*!
             *   piston theory boundary condition for the whole domain
             */
            MAST::TimeDomainFlutterSolver*                        _flutter_solver;
            
            /*!
             *   flutter root from the analysis
             */
            MAST::FlutterRootBase*                               _flutter_root;

            // vector of basis vectors from modal analysis
            std::vector<libMesh::NumericVector<Real>*>           _basis;
        };
    }
}

#endif // __mast__structural_example_base_h__
