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

#ifndef __mast_plate_thermal_stress_optimization_h__
#define __mast_plate_thermal_stress_optimization_h__

// C++ includes
#include <memory>

// MAST includes
#include "examples/base/multilinear_interpolation.h"
#include "examples/structural/plate_optimization/plate_optimization_base.h"
#include "base/field_function_base.h"
#include "base/physics_discipline_base.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/structural_system_initialization.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "optimization/gcmma_optimization_interface.h"
#include "optimization/function_evaluation.h"
#include "boundary_condition/dirichlet_boundary_condition.h"


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/serial_mesh.h"
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
    class StructuralSystemInitialization;
    class StructuralDiscipline;
    class Parameter;
    class ConstantFieldFunction;
    class IsotropicMaterialPropertyCard;
    class Solid2DSectionElementPropertyCard;
    class DirichletBoundaryCondition;
    class BoundaryConditionBase;
    class StressStrainOutputBase;
    class StructuralNonlinearAssembly;
    class SectionOffset;
    
    
    struct PlateBendingThermalStressSizingOptimization:
    public MAST::FunctionEvaluation {
        
        
        PlateBendingThermalStressSizingOptimization
        (const libMesh::Parallel::Communicator& comm);
        
        
        ~PlateBendingThermalStressSizingOptimization();
        
        /*!
         *   initializes the object for specified characteristics
         */
        void init(GetPot& infile, libMesh::ElemType e_type, bool if_vk);

        
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
         *  @returns a pointer to the function that evaluates the objective
         */
        virtual MAST::FunctionEvaluation::funobj
        get_objective_evaluation_function();
        
        
        /*!
         *  @returns a pointer to the function that evaluates the constraint
         */
        virtual MAST::FunctionEvaluation::funcon
        get_constraint_evaluation_function();

        /*!
         *   clears the stress data structures for a followup analysis
         */
        void clear_stresss();
        
        
        bool _initialized;

        
        // length of domain
        Real _length;
        
        // width of domain
        Real _width;

        // length of domain
        Real _stress_limit;
        
        // number of elements and number of stations at which DVs are defined
        unsigned int
        _n_divs_x,
        _n_divs_y,
        _n_elems,
        _n_stations_x;
        
        // create the mesh
        libMesh::SerialMesh*           _mesh;
        
        // create the equation system
        libMesh::EquationSystems*      _eq_sys;
        
        // create the libmesh system
        libMesh::NonlinearImplicitSystem*  _sys;
        
        // initialize the system to the right set of variables
        MAST::StructuralSystemInitialization* _structural_sys;
        MAST::StructuralDiscipline*           _discipline;
        
        // nonlinear assembly object
        MAST::StructuralNonlinearAssembly *_assembly;
        
        // create the property functions and add them to the
        MAST::Parameter
        *_E,
        *_alpha,
        *_nu,
        *_kappa,
        *_rho,
        *_temp,
        *_zero;
        
        MAST::ConstantFieldFunction
        *_E_f,
        *_nu_f,
        *_alpha_f,
        *_kappa_f,
        *_rho_f,
        *_temp_f,
        *_ref_temp_f;
        
        // Section offset
        MAST::SectionOffset*                       _hoff_f;
        
        // Weight function to calculate the weight of the structure
        MAST::PlateWeight *_weight;
        
        // create the material property card
        MAST::IsotropicMaterialPropertyCard*            _m_card;
        
        // create the element property card
        MAST::Solid2DSectionElementPropertyCard*        _p_card;
        
        // create the Dirichlet boundary condition on left edge
        MAST::DirichletBoundaryCondition*               _dirichlet_left;
        
        // create the Dirichlet boundary condition on right edge
        MAST::DirichletBoundaryCondition*               _dirichlet_right;
        
        // create the Dirichlet boundary condition on bottom edge
        MAST::DirichletBoundaryCondition*               _dirichlet_bottom;
        
        // create the Dirichlet boundary condition on top edge
        MAST::DirichletBoundaryCondition*               _dirichlet_top;
        
        // create the temperature load
        MAST::BoundaryConditionBase*             _T_load;
        
        // output quantity objects to evaluate stress
        std::vector<MAST::StressStrainOutputBase*>      _outputs;
        
        
        // stationwise parameter definitions
        std::vector<MAST::Parameter*>                   _th_station_parameters;
        
        // stationwise function objects for thickness
        std::vector<MAST::ConstantFieldFunction*>       _th_station_functions;
        
        /*!
         *   interpolates thickness between stations
         */
        std::auto_ptr<MAST::MultilinearInterpolation>   _th_f;
        
        /*!
         *   scaling parameters for design optimization problem
         */
        std::vector<Real>
        _dv_scaling,
        _dv_low,
        _dv_init;
    };
}


#endif /* __mast_plate_thermal_stress_optimization_h__ */

