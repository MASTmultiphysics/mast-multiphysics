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

#ifndef __mast_beam_euler_fsi_h__
#define __mast_beam_euler_fsi_h__


// C++ includes
#include <memory>
#include <vector>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/fe_type.h"
#include "libmesh/dof_map.h"
#include "libmesh/parallel.h"



namespace MAST {
    
    // Forward declerations
    class ConservativeFluidSystemInitialization;
    class ConservativeFluidDiscipline;
    class StructuralSystemInitialization;
    class Parameter;
    class PhysicsDisciplineBase;
    class ConstantFieldFunction;
    class IsotropicMaterialPropertyCard;
    class Solid1DSectionElementPropertyCard;
    class DirichletBoundaryCondition;
    class BoundaryConditionBase;
    class NonlinearSystem;
    class FlightCondition;
    class FrequencyFunction;
    template <typename ValType> class FieldFunction;
    class MeshFieldFunction;
    class NormalRotationMeshFunction;
    class PressureFunction;
    class FrequencyDomainPressureFunction;
    class MultiphysicsNonlinearSolverBase;
    class FSIBoundaryConditionUpdates;
    class AugmentGhostElementSendListObj;
    
    struct BeamEulerFSIAnalysis {
        
        
        BeamEulerFSIAnalysis();
        
        
        ~BeamEulerFSIAnalysis();
        
        
        /*!
         *   @returns a pointer to the parameter of the specified name.
         *   If no parameter exists by the specified name, then a \p nullptr
         *   pointer is returned and a message is printed with a valid list
         *   of parameters.
         */
        MAST::Parameter* get_parameter(const std::string& nm);
        
        
        /*!
         *  solves the system and returns the flutter velocity
         */
        Real solve(bool if_write_output = false,
                   const Real tol = 1.e-1,
                   const unsigned int max_bisection_iters = 20);
        
        
        /*!
         *  solves the sensitivity of system and returns the sensitiivty of
         *  flutter speed
         */
        Real sensitivity_solve(const MAST::Parameter& p);
        
        // parameters to control the time step size
        unsigned int _max_time_steps;
        Real         _time_step_size;
        
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
        MAST::PhysicsDisciplineBase*              _structural_discipline;
        
        
        // initialize the system to the right set of variables
        MAST::ConservativeFluidSystemInitialization* _fluid_sys_init;
        MAST::ConservativeFluidDiscipline*           _fluid_discipline;
        
        
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
        MAST::MeshFieldFunction                    *_vel;
        MAST::MeshFieldFunction                    *_displ;
        MAST::NormalRotationMeshFunction           *_normal_rot;
        
        
        /*!
         *   surface pressure
         */
        MAST::PressureFunction                *_pressure_function;
        
        
        // parameters used in the system
        MAST::Parameter
        *_velocity;
        
        
        MAST::ConstantFieldFunction
        *_velocity_f;
        
        
        Real
        _length;
        
        
        // create the property functions and add them to the
        MAST::Parameter
        *_thy,
        *_thz,
        *_rho,
        *_E,
        *_nu,
        *_alpha,
        *_temp,
        *_zero;
        
        MAST::ConstantFieldFunction
        *_thy_f,
        *_thz_f,
        *_rho_f,
        *_E_f,
        *_nu_f,
        *_alpha_f,
        *_temp_f,
        *_ref_temp_f,
        *_hyoff_f,
        *_hzoff_f;
        
        
        // create the material property card
        MAST::IsotropicMaterialPropertyCard*     _m_card;
        
        // create the element property card
        MAST::Solid1DSectionElementPropertyCard* _p_card;
        
        // create the Dirichlet boundary condition on left edge
        MAST::DirichletBoundaryCondition*        _dirichlet_left;
        
        // create the Dirichlet boundary condition on right edge
        MAST::DirichletBoundaryCondition*        _dirichlet_right;
        
        // create the temperature load
        MAST::BoundaryConditionBase*             _T_load;

        // object to augment the send list of ghosted fluid elements
        MAST::AugmentGhostElementSendListObj*    _augment_send_list_obj;

        // updates the boundary condition objects with the current solution
        MAST::FSIBoundaryConditionUpdates*      _bc_updates;
        
        // vector of parameters to evaluate sensitivity wrt
        std::vector<MAST::Parameter*>           _params_for_sensitivity;
    };
}



#endif //  __mast_beam_euler_fsi_h__

