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

#ifndef __mast_stiffened_plate_thermally_stressed_flutter_optimization_h__
#define __mast_stiffened_plate_thermally_stressed_flutter_optimization_h__

// C++ includes
#include <memory>

// MAST includes
#include "examples/base/multilinear_interpolation.h"
#include "base/field_function_base.h"
#include "base/physics_discipline_base.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/structural_system_initialization.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "optimization/gcmma_optimization_interface.h"
#include "optimization/function_evaluation.h"
#include "boundary_condition/dirichlet_boundary_condition.h"


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/fe_type.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh_function.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/getpot.h"


// get this from the global namespace
extern libMesh::LibMeshInit* _init;


namespace MAST {
    
    
    // Forward declerations
    class ConservativeFluidSystemInitialization;
    class ConservativeFluidDiscipline;
    class StructuralSystemInitialization;
    class StructuralDiscipline;
    class Parameter;
    class ConstantFieldFunction;
    class IsotropicMaterialPropertyCard;
    class Solid1DSectionElementPropertyCard;
    class Solid2DSectionElementPropertyCard;
    class DirichletBoundaryCondition;
    class BoundaryConditionBase;
    class StructuralNonlinearAssembly;
    class UGFlutterSolver;
    class FlutterRootBase;
    class FlightCondition;
    class FrequencyFunction;
    class ComplexMeshFieldFunction;
    class ComplexNormalRotationMeshFunction;
    class PressureFunction;
    class FrequencyDomainPressureFunction;
    class NonlinearSystem;
    class StructuralModalEigenproblemAssembly;
    class FSIGeneralizedAeroForceAssembly;
    class FrequencyDomainLinearizedComplexAssembly;
    class ComplexSolverBase;
    class GAFDatabase;
    class AugmentGhostElementSendListObj;
    class StiffenedPlateWeight;
    class StructuralNearNullVectorSpace;
    class StressStrainOutputBase;
    
    struct StiffenedPlateThermallyStressedFSIFlutterSizingOptimization:
    public MAST::FunctionEvaluation {
        
        
        StiffenedPlateThermallyStressedFSIFlutterSizingOptimization
        (const libMesh::Parallel::Communicator& comm);
        
        
        ~StiffenedPlateThermallyStressedFSIFlutterSizingOptimization();
        
        
        void init(GetPot &infile, libMesh::ElemType etype, bool if_nonlin);

        
        /*!
         *   initialize the design variables values and bounds
         */
        virtual void init_dvar(std::vector<Real>& x,
                               std::vector<Real>& xmin,
                               std::vector<Real>& xmax);
        
        
        /*!
         *    the core routine that performs the function evaluations
         */
        virtual void evaluate(const std::vector<Real>& dvars,
                              Real& obj,
                              bool eval_obj_grad,
                              std::vector<Real>& obj_grad,
                              std::vector<Real>& fvals,
                              std::vector<bool>& eval_grads,
                              std::vector<Real>& grads);
        
        /*!
         *   customized output
         */
        virtual void output(unsigned int iter,
                            const std::vector<Real>& x,
                            Real obj,
                            const std::vector<Real>& fval,
                            bool if_write_to_optim_file) const;

        /*!
         *   clears the stored stress data, which is done before each 
         *   stress evaluation
         */
        void clear_stresss();
        
        
        /*!
         *  @returns a pointer to the function that evaluates the objective
         */
        virtual MAST::FunctionEvaluation::funobj
        get_objective_evaluation_function();
        
        
        /*!
         *  @returns a pointer to the function that evaluates the constraint
         */
        virtual MAST::FunctionEvaluation::funcon
        get_constraint_evaluation_function();
        
        
        bool _initialized;
        
        Real
        _length,
        _width,
        _stress_limit,
        _k_lower,
        _k_upper,
        _V0_flutter;

        // number of elements and number of stations at which DVs are defined
        unsigned int
        _n_eig,
        _n_elems,
        _n_divs_x,
        _n_divs_between_stiff,
        _n_stiff, // assumed along x-axis with equidistance spacing
        _n_plate_elems,
        _n_elems_per_stiff,
        _n_stations,
        _n_k_divs,
        _n_dv_stations_x,
        _n_load_steps;

        
        // create the structural mesh
        libMesh::SerialMesh*                     _structural_mesh;
        
        
        // create the fluid mesh
        libMesh::ParallelMesh*                     _fluid_mesh;
        
        
        // create the equation system
        libMesh::EquationSystems
        *_structural_eq_sys,
        *_fluid_eq_sys;
        
        
        // create the libmesh system
        MAST::NonlinearSystem
        *_structural_sys,
        *_fluid_sys;

        // initialize the system to the right set of variables
        MAST::StructuralSystemInitialization*    _structural_sys_init;
        MAST::StructuralDiscipline*              _structural_discipline;
        
        MAST::StructuralNearNullVectorSpace*       _nsp;
        
        // initialize the system to the right set of variables
        MAST::ConservativeFluidSystemInitialization* _fluid_sys_init;
        MAST::ConservativeFluidDiscipline*           _fluid_discipline;
        
        
        // frequency domain assembly and solver objects
        MAST::FrequencyDomainLinearizedComplexAssembly   *_frequency_domain_fluid_assembly;
    
        MAST::ComplexSolverBase                          *_complex_solver;
        
        // modal assembly object
        MAST::StructuralModalEigenproblemAssembly    *_modal_assembly;
        
        // structural nonlinear assembly object
        MAST::StructuralNonlinearAssembly     *_structural_nonlinear_assembly;

        // flight condition
        MAST::FlightCondition*                       _flight_cond;
        
        
        // boundary condition
        MAST::BoundaryConditionBase
        *_far_field,
        *_symm_wall,
        *_slip_wall,
        *_pressure;
        
        
        /*!
         *   surface motion
         */
        MAST::ComplexMeshFieldFunction                    *_displ;
        MAST::ComplexNormalRotationMeshFunction           *_normal_rot;
        
        /*!
         *   surface pressure
         */
        MAST::PressureFunction                *_pressure_function;
        MAST::FrequencyDomainPressureFunction *_freq_domain_pressure_function;
        
        // parameters used in the system
        MAST::Parameter
        *_omega,
        *_velocity,
        *_b_ref;
        
        
        MAST::ConstantFieldFunction
        *_omega_f,
        *_velocity_f,
        *_b_ref_f;
        
        
        /*!
         *   frequency object
         */
        MAST::FrequencyFunction
        *_freq_function;

        /*!
         *   flutter solver
         */
        MAST::UGFlutterSolver*                   _flutter_solver;
        
        // vector of basis vectors from modal analysis
        std::vector<libMesh::NumericVector<Real>*>
        _aero_basis,
        _structural_basis;

        
        // map of reduced order matrices corresponding to the basis, stored
        // for multiple reduced frequency values. Linear interpolation
        // is used for evaluation
        MAST::GAFDatabase*   _gaf_database;
        
        // create the property functions and add them to the
        MAST::Parameter
        *_E,
        *_nu,
        *_kappa,
        *_alpha,
        *_rho,
        *_temp,
        *_p_cav,
        *_zero;
        
        MAST::ConstantFieldFunction
        *_E_f,
        *_nu_f,
        *_kappa_f,
        *_alpha_f,
        *_rho_f,
        *_hoff_plate_f,
        *_thzoff_stiff_f,
        *_temp_f,
        *_ref_temp_f,
        *_p_cav_f;
        

        // section offset for each stiffener
        std::vector<MAST::SectionOffset*>
        _hyoff_stiff_f,
        _hzoff_stiff_f;

        // Weight function to calculate the weight of the structure
        MAST::StiffenedPlateWeight*                     _weight;
        
        // create the material property card
        MAST::IsotropicMaterialPropertyCard*            _m_card;
        
        // create the element property card
        MAST::Solid2DSectionElementPropertyCard*        _p_card_plate;

        // create the element property card, one for each stiffener
        std::vector<MAST::Solid1DSectionElementPropertyCard*>   _p_card_stiff;

        // create the Dirichlet boundary condition
        MAST::DirichletBoundaryCondition
        *_dirichlet_left,
        *_dirichlet_right,
        *_dirichlet_bottom,
        *_dirichlet_top;

        MAST::BoundaryConditionBase
        *_T_load,
        *_p_load;

        // output quantity objects to evaluate stress
        std::vector<MAST::StressStrainOutputBase*>     _outputs;

        // object to augment the send list of ghosted fluid elements
        MAST::AugmentGhostElementSendListObj*    _augment_send_list_obj;

        // stationwise parameter definitions
        std::vector<MAST::Parameter*>
        _th_station_parameters_plate,
        _thy_station_parameters_stiff,
        _thz_station_parameters_stiff; // N_stiff*N_stations params

        
        // stationwise function objects for thickness
        std::vector<MAST::ConstantFieldFunction*>
        _th_station_functions_plate,
        _thy_station_functions_stiff,
        _thz_station_functions_stiff; // N_stiff*N_stations functions;
        
        /*!
         *   interpolates thickness between stations
         */
        std::unique_ptr<MAST::MultilinearInterpolation>   _th_plate_f;
        
        /*!
         *   interpolates thickness between stations
         */
        std::vector<MAST::MultilinearInterpolation*>
        _thy_stiff_f,
        _thz_stiff_f; // one per stiffener

        /*!
         *   vector of problem parameters
         */
        std::vector<MAST::Parameter*>                _problem_parameters;

        /*!
         *   scaling parameters for design optimization problem
         */
        std::vector<Real>
        _dv_scaling,
        _dv_low,
        _dv_init;
    };
}


#endif // __mast_stiffened_plate_thermally_stressed_flutter_optimization_h__
