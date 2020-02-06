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

#ifndef __mast_thermoelasticity_topology_optimization_level_set_2d_h__
#define __mast_thermoelasticity_topology_optimization_level_set_2d_h__


// MAST includes
#include "examples/structural/topology_optim_2D/topology_optim_2D.h"


namespace MAST  {

    namespace Examples {

        // Forward declerations
        class TemperatureFunction;
        
        class ThermoelasticityTopologyOptimizationLevelSet2D:
        public MAST::Examples::TopologyOptimizationLevelSet2D {

        public:

            ThermoelasticityTopologyOptimizationLevelSet2D(const libMesh::Parallel::Communicator& comm_in);

            virtual ~ThermoelasticityTopologyOptimizationLevelSet2D();

            /*!
             *   \p grads(k): Derivative of f_i(x) with respect
             *   to x_j, where k = (j-1)*M + i.
             */
            virtual void evaluate(const std::vector<Real>& dvars,
                                  Real& obj,
                                  bool eval_obj_grad,
                                  std::vector<Real>& obj_grad,
                                  std::vector<Real>& fvals,
                                  std::vector<bool>& eval_grads,
                                  std::vector<Real>& grads);

        protected:

            virtual void _init_system_and_discipline();
            virtual void _init_dirichlet_conditions();
            virtual void _init_loads();
            virtual void _init_temperature_load();
            virtual void _init_section_property();
            virtual void _init_conduction_boundary_dirichlet_constraint(const unsigned int bid,
                                                                        const std::string& tag);
            void _evaluate_constraint_sensitivity
            (MAST::StressStrainOutputBase& stress,
             MAST::AssemblyElemOperations& conduction_elem_ops,
             MAST::AssemblyElemOperations& nonlinear_elem_ops,
             MAST::LevelSetNonlinearImplicitAssembly& conduction_assembly,
             MAST::LevelSetNonlinearImplicitAssembly& nonlinear_assembly,
             MAST::StructuralModalEigenproblemAssemblyElemOperations& eigen_elem_ops,
             MAST::LevelSetEigenproblemAssembly& eigen_assembly,
             const std::vector<bool>& eval_grads,
             std::vector<Real>& grads);
            
            MAST::NonlinearSystem*                    _conduction_sys;
            MAST::HeatConductionSystemInitialization* _conduction_sys_init;
            MAST::PhysicsDisciplineBase*              _conduction_discipline;
            MAST::Examples::TemperatureFunction*      _temp_function;
        };
    }
}


#endif //  __mast_thermoelasticity_topology_optimization_level_set_2d_h__


