/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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

#ifndef __mast_beam_optimization_h__
#define __mast_beam_optimization_h__

// C++ includes
#include <memory>

// MAST includes
#include "base/field_function_base.h"
#include "base/physics_discipline_base.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/structural_system_initialization.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_1d_section_element_property_card.h"
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
    
    /*!
     *   This class provides the ability to interpolate a function in between
     *   a set of tabulated points.
     */
    class MultilinearInterpolation:
    public MAST::FieldFunction<Real> {
    public:
        MultilinearInterpolation(const std::string& nm,
                                 std::map<Real, MAST::FieldFunction<Real>*>& values);
        
        MultilinearInterpolation(const MAST::MultilinearInterpolation& o);
        
        
        virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const;
        
        
        virtual ~MultilinearInterpolation();
        
    protected:
        
        std::map<Real, MAST::FieldFunction<Real>*> _values;
        
    public:
        
        virtual void operator() (const libMesh::Point& p, Real t, Real& v) const;
        
        virtual void derivative(const MAST::DerivativeType d,
                                const MAST::FunctionBase& f,
                                const libMesh::Point& p,
                                Real t,
                                Real& v) const;
    };
    
    
    /*!
     *   Function object evaluates the beam offset for the specified height
     */
    class BeamOffset: public MAST::FieldFunction<Real> {
    public:
        BeamOffset(const std::string& nm,
                   MAST::FieldFunction<Real> *thickness);
        
        BeamOffset(const MAST::BeamOffset& o);
        
        virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const;
        
        virtual ~BeamOffset();
        
    protected:
        
        MAST::FieldFunction<Real> *_dim;
        
    public:
        
        virtual void operator() (const libMesh::Point& p, Real t, Real& v) const;
        
        virtual void derivative(const MAST::DerivativeType d,
                                const MAST::FunctionBase& f,
                                const libMesh::Point& p,
                                Real t,
                                Real& v) const;
        
    };
    
    
    
    /*!
     *   Function object evaluates the weight and its sensitivity with
     *   respect to the specified variable.
     */
    class Weight: public MAST::FieldFunction<Real> {
    public:
        
        /*!
         *   Constructor requires the mesh and the
         */
        Weight(MAST::PhysicsDisciplineBase& discipline);
        
        /*!
         *  copy constructor
         */
        Weight(const MAST::Weight& w);
        
        /*!
         *  @returns a new object as a clone, encapsulated in a smart-pointer
         */
        virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const;
        
        /*!
         *  virtual destructor
         */
        virtual ~Weight();
        
    protected:
        
        /*!
         *  discipline object provides the mesh and material properties for
         *  calculation of the mass
         */
        MAST::PhysicsDisciplineBase& _discipline;
        
    public:
        
        /*!
         *   overloaded operator evaluates and returns the mass of the given
         *   structural model.
         */
        virtual void operator() (const libMesh::Point& p,
                                 Real t,
                                 Real& v) const ;
        
        
        
        /*!
         *   evaluates the sensitivity of structural mass with respect
         *   to the design variable.
         */
        virtual void derivative(const MAST::DerivativeType d,
                                const MAST::FunctionBase& f,
                                const libMesh::Point& p,
                                Real t,
                                Real& v) const ;
    };
    
    
    // Forward declerations
    class StructuralSystemInitialization;
    class StructuralDiscipline;
    class Parameter;
    class ConstantFieldFunction;
    class IsotropicMaterialPropertyCard;
    class Solid1DSectionElementPropertyCard;
    class DirichletBoundaryCondition;
    class BoundaryConditionBase;
    class StressStrainOutputBase;
    class StructuralNonlinearAssembly;
    
    
    struct BeamBendingSizingOptimization:
    public MAST::FunctionEvaluation {
        
        
        BeamBendingSizingOptimization(GetPot& infile,
                                      std::ostream& output);
        
        
        ~BeamBendingSizingOptimization();
        
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
        
        
        /*!
         *   clears the stress data structures for a followup analysis
         */
        void clear_stresss();
        
        
        // length of domain
        Real _length;

        
        // length of domain
        Real _stress_limit;

        // number of elements and number of stations at which DVs are defined
        unsigned int
        _n_elems,
        _n_stations;
        
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
        *_thz,
        *_E,
        *_nu,
        *_rho,
        *_press,
        *_zero;
        
        MAST::ConstantFieldFunction
        *_thz_f,
        *_E_f,
        *_nu_f,
        *_rho_f,
        *_hyoff_f,
        *_hzoff_f,
        *_press_f;
        
        
        // Weight function to calculate the weight of the structure
        MAST::Weight *_weight;
        
        // create the material property card
        MAST::IsotropicMaterialPropertyCard*            _m_card;
        
        // create the element property card
        MAST::Solid1DSectionElementPropertyCard*        _p_card;
        
        // create the Dirichlet boundary condition on left edge
        MAST::DirichletBoundaryCondition*               _dirichlet_left;
        
        // create the Dirichlet boundary condition on right edge
        MAST::DirichletBoundaryCondition*               _dirichlet_right;
        
        // create the pressure boundary condition
        MAST::BoundaryConditionBase*                    _p_load;
        
        // output quantity objects to evaluate stress
        std::vector<MAST::StressStrainOutputBase*>      _outputs;
        

        // stationwise parameter definitions
        std::vector<MAST::Parameter*>                   _thy_station_parameters;
        
        // stationwise function objects for thickness
        std::vector<MAST::ConstantFieldFunction*>       _thy_station_functions;
        
        /*!
         *   interpolates thickness between stations
         */
        std::auto_ptr<MAST::MultilinearInterpolation>   _thy_f;
        
        /*!
         *   scaling parameters for design optimization problem
         */
        std::vector<Real>
        _dv_scaling,
        _dv_low,
        _dv_init;
    };
}


#endif /* __mast_beam_optimization_h__ */
