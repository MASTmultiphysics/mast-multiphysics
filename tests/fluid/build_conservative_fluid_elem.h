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


#ifndef __mast_build_conservative_fluid_elem_h__
#define __mast_build_conservative_fluid_elem_h__


// C++ includes
#include <memory>
#include <vector>
#include <string>

// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/fe_type.h"
#include "libmesh/dof_map.h"



namespace MAST {
    
    // Forward declerations
    class ConservativeFluidSystemInitialization;
    class ConservativeFluidDiscipline;
    class Parameter;
    class ConstantFieldFunction;
    class BoundaryConditionBase;
    class FlightCondition;
    class FrequencyFunction;
    class RigidSurfaceMotion;
    class RigidSurfaceDisplacement;
    class RigidSurfaceNormalRotation;
    class NonlinearSystem;
    
    
    struct BuildConservativeFluidElem {
        
        
        BuildConservativeFluidElem();
        
        
        ~BuildConservativeFluidElem();
        
        
        // create the mesh
        libMesh::ParallelMesh*           _mesh;
        
        // create the equation system
        libMesh::EquationSystems*      _eq_sys;
        
        // create the libmesh system
        MAST::NonlinearSystem*  _sys;
        
        // initialize the system to the right set of variables
        MAST::ConservativeFluidSystemInitialization* _fluid_sys;
        MAST::ConservativeFluidDiscipline*           _discipline;
        
        // flight condition
        MAST::FlightCondition*          _flight_cond;
        
        
        RealVectorX _base_sol;
        
        // boundary condition
        MAST::BoundaryConditionBase
        *_far_field,
        *_slip_wall;
        
        /*!
         *   surface rigid motion
         */
        MAST::RigidSurfaceMotion
        *_motion;
        MAST::RigidSurfaceDisplacement
        *_displacement;
        MAST::RigidSurfaceNormalRotation
        *_normal_rot;
        
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
        MAST::FrequencyFunction        *_freq_function;
    };
}




#endif // __mast_build_conservative_fluid_elem_h__

