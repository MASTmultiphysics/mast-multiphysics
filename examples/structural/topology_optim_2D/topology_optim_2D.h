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

#ifndef __mast_topology_optimization_level_set_2d_h__
#define __mast_topology_optimization_level_set_2d_h__


// MAST includes
#include "examples/structural/base/structural_example_2d.h"
#include "optimization/function_evaluation.h"


namespace MAST  {
    
    // Forward declerations
    class LevelSetSystemInitialization;
    class LevelSetDiscipline;
    class LevelSetVolume;
    class LevelSetNonlinearImplicitAssembly;
    class LevelSetEigenproblemAssembly;
    class StressStrainOutputBase;
    class LevelSetBoundaryVelocity;
    class PhiMeshFunction;
    class AssemblyElemOperations;
    class StructuralModalEigenproblemAssemblyElemOperations;
    class HeatConductionSystemInitialization;
    template <typename ValType> class FieldFunction;
    
    
    namespace Examples {
        
        class TopologyOptimizationLevelSet2D:
        public MAST::Examples::StructuralExample2D,
        public MAST::FunctionEvaluation {
            
        public:
            
            TopologyOptimizationLevelSet2D(const libMesh::Parallel::Communicator& comm_in);
            
            virtual ~TopologyOptimizationLevelSet2D();

            virtual void initialize_solution();
            
            /*!
             *    initializes the design data after calling the parent
             *    class' init method
             */
            virtual void init(MAST::Examples::GetPotWrapper& input,
                              const std::string& prefix);

            
            /*!
             *   initializes the design variable vector, called by the
             *   optimization interface.
             */
            virtual void init_dvar(std::vector<Real>& x,
                                   std::vector<Real>& xmin,
                                   std::vector<Real>& xmax);
            
            /*!
             *   \par grads(k): Derivative of f_i(x) with respect
             *   to x_j, where k = (j-1)*M + i.
             */
            virtual void evaluate(const std::vector<Real>& dvars,
                                  Real& obj,
                                  bool eval_obj_grad,
                                  std::vector<Real>& obj_grad,
                                  std::vector<Real>& fvals,
                                  std::vector<bool>& eval_grads,
                                  std::vector<Real>& grads);

            /*!
             *   solves the level set equations for either propagation or
             *   reinitialization of the level set function
             */
            void level_set_solve();
            

            virtual void output(unsigned int iter,
                                const std::vector<Real>& x,
                                Real obj,
                                const std::vector<Real>& fval,
                                bool if_write_to_optim_file);

        protected:
            
            virtual void _init_mesh();
            virtual void _init_system_and_discipline();
            virtual void _init_dirichlet_conditions();
            virtual void _init_eq_sys();
            virtual void _init_loads();
            virtual void _init_material();
            virtual void _init_section_property();
            virtual void _init_phi_dvs();
            void _evaluate_volume_sensitivity(LevelSetVolume& volume,
                                              MAST::LevelSetNonlinearImplicitAssembly& assembly,
                                              std::vector<Real>& obj_grad);
            void _evaluate_constraint_sensitivity
            (MAST::StressStrainOutputBase& stress,
             MAST::AssemblyElemOperations& nonlinear_elem_ops,
             MAST::LevelSetNonlinearImplicitAssembly& nonlinear_assembly,
             MAST::StructuralModalEigenproblemAssemblyElemOperations& eigen_elem_ops,
             MAST::LevelSetEigenproblemAssembly& eigen_assembly,
             const std::vector<bool>& eval_grads,
             std::vector<Real>& grads);
            
            Real                                      _stress_lim;
            Real                                      _p_val;
            Real                                      _ref_eig_val;
            unsigned int                              _n_eig_vals;
            libMesh::FEType                           _level_set_fetype;
            libMesh::UnstructuredMesh*                _level_set_mesh;
            libMesh::EquationSystems*                 _level_set_eq_sys;
            MAST::NonlinearSystem*                    _level_set_sys;
            MAST::NonlinearSystem*                    _level_set_sys_on_str_mesh;
            MAST::NonlinearSystem*                    _indicator_sys;
            MAST::LevelSetSystemInitialization*       _level_set_sys_init_on_str_mesh;
            MAST::LevelSetSystemInitialization*       _level_set_sys_init;
            MAST::HeatConductionSystemInitialization* _indicator_sys_init;
            MAST::PhysicsDisciplineBase*              _indicator_discipline;
            MAST::LevelSetDiscipline*                 _level_set_discipline;
            PhiMeshFunction*                          _level_set_function;
            MAST::LevelSetBoundaryVelocity*           _level_set_vel;
            libMesh::ExodusII_IO*                     _output;
            std::vector<std::pair<unsigned int, MAST::Parameter*>>  _dv_params;
        };
    }
}


#endif //  __mast_topology_optimization_level_set_2d_h__

