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
#include "base/field_function_base.h"
#include "base/mesh_field_function.h"


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
    
        
        class FluxLoad:
        public MAST::FieldFunction<Real> {
        public:
            FluxLoad(const std::string& nm, Real p, Real l1, Real fraction):
            MAST::FieldFunction<Real>(nm), _p(p), _l1(l1), _frac(fraction) { }
            virtual ~FluxLoad() {}
            virtual void operator() (const libMesh::Point& p, const Real t, Real& v) const {
                if (fabs(p(0)-_l1*0.5) <= 0.5*_frac*_l1) v = _p;
                else v = 0.;
            }
            virtual void derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, Real& v) const {
                v = 0.;
            }
        protected:
            Real _p, _l1, _frac;
        };
        

        
        class PhiMeshFunction:
        public MAST::FieldFunction<Real> {
        public:
            PhiMeshFunction():
            MAST::FieldFunction<Real>("phi"), _phi(nullptr) { }
            virtual ~PhiMeshFunction(){ if (_phi) delete _phi;}
            
            void init(MAST::SystemInitialization& sys, const libMesh::NumericVector<Real>& sol) {
                if (!_phi) _phi = new MAST::MeshFieldFunction(sys, "phi");
                else _phi->clear();
                _phi->init(sol);
            }
            
            MAST::MeshFieldFunction& get_mesh_function() {return *_phi;}
            
            virtual void operator() (const libMesh::Point& p, const Real t, Real& v) const {
                libmesh_assert(_phi);
                RealVectorX v1;
                (*_phi)(p, t, v1);
                v = v1(0);
            }
            
        protected:
            MAST::MeshFieldFunction *_phi;
        };

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
            virtual void _init_indicator_system_dirichlet_conditions();
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
            Real                                      _p_val, _vm_rho;
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

