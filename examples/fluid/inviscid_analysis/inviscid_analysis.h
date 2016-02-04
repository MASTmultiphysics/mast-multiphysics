/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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

#ifndef __mast_inviscid_analysis_h__
#define __mast_inviscid_analysis_h__


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
#include "libmesh/nonlinear_implicit_system.h"
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
    
    
    struct InviscidAnalysis {
        
        
        InviscidAnalysis();
        
        
        ~InviscidAnalysis();
        
        
        /*!
         *   @returns a pointer to the parameter of the specified name.
         *   If no parameter exists by the specified name, then a \p NULL
         *   pointer is returned and a message is printed with a valid list
         *   of parameters.
         */
        MAST::Parameter* get_parameter(const std::string& nm);
        
        /*!
         *  solves the system and returns the final solution
         */
        const libMesh::NumericVector<Real>&
        solve(bool if_write_output = false);
        
        
        /*!
         *  solves the sensitivity of system and returns the final solution
         */
        const libMesh::NumericVector<Real>&
        sensitivity_solve(MAST::Parameter& p,
                          bool if_write_output = false);
        
        // parameters to control the time step size
        unsigned int _max_time_steps;
        Real         _time_step_size;
        
        
        // create the mesh
        libMesh::ParallelMesh*           _mesh;
        
        // create the equation system
        libMesh::EquationSystems*      _eq_sys;
        
        // create the libmesh system
        libMesh::NonlinearImplicitSystem*  _sys;
        
        // initialize the system to the right set of variables
        MAST::ConservativeFluidSystemInitialization* _fluid_sys;
        MAST::ConservativeFluidDiscipline*           _discipline;

        // flight condition
        MAST::FlightCondition*          _flight_cond;

        MAST::BoundaryConditionBase
        *_far_field,
        *_slip_wall;

        
        // vector of parameters to evaluate sensitivity wrt
        std::vector<MAST::Parameter*> _params_for_sensitivity;
    };
}



#endif //  __mast_inviscid_analysis_h__
