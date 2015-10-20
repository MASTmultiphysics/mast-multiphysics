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
#include "boundary_condition/dirichlet_boundary_condition.h"


// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/mesh_function.h"
#include "libmesh/parameter_vector.h"


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
    
    
    /*!
     *   This class implements the function evaluation routine that provides
     *   values and gradients of the objective and constraint functions
     *   for this optimization problem
     */
    class SizingOptimization:
    public MAST::FunctionEvaluation {
        
    public:
        
        SizingOptimization(libMesh::LibMeshInit& init,
                           std::ostream& output);
        
        virtual ~SizingOptimization();
        
        
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
        
    protected:
        
        /*!
         *   member function that initializes the data structures. Gets called
         *   by the constructor.
         */
        void _init();
        
        
        
        
        libMesh::LibMeshInit                            &_libmesh_init;
        unsigned int
        _n_elems,
        _n_stations;
        
        Real
        _disp_0,
        _vf_0;
        
        MAST::Weight                                    *_weight;
        
        libMesh::UnstructuredMesh                       *_mesh;
        libMesh::EquationSystems                        *_eq_systems;
        libMesh::NonlinearImplicitSystem                *_static_system;
        std::auto_ptr<MAST::StructuralDiscipline>       _structural_discipline;
        std::auto_ptr<MAST::StructuralSystemInitialization> _structural_sys;
        std::auto_ptr<MAST::DirichletBoundaryCondition> _dirichlet;
        std::auto_ptr<MAST::BoundaryConditionBase>      _flux_load;
        std::auto_ptr<MAST::MaterialPropertyCardBase>   _materials;
        std::auto_ptr<MAST::ElementPropertyCardBase>    _elem_properties;
        std::auto_ptr<libMesh::MeshFunction>            _disp_function;
        std::vector<libMesh::MeshFunction*>             _disp_function_sens;
        std::auto_ptr<MAST::Parameter>
        _pressure,
        _E,
        _nu,
        _rho,
        _kappa,
        _h_z,
        _offset_h_z,
        _temperature,
        _ref_temperature;
        std::auto_ptr<MAST::ConstantFieldFunction>
        _pressure_fn,
        _E_fn,
        _nu_fn,
        _rho_fn,
        _kappa_fn,
        _h_z_fn,
        _offset_h_z_fn,
        _temperature_fn,
        _ref_temperature_fn;
        
        /*!
         *    map of thickness functions at discrete stations
         */
        std::map<Real, MAST::FieldFunction<Real>*>      _h_y_station_vals;
        std::auto_ptr<MAST::MultilinearInterpolation>   _h_y_fn;
        std::auto_ptr<MAST::BeamOffset>                 _offset_h_y_fn;
        
        
        libMesh::ParameterVector                        _parameters;
        std::vector<MAST::ConstantFieldFunction*>       _parameter_functions;
        
        
        std::vector<Real>
        _dv_scaling,
        _dv_low,
        _dv_init;
    };
}


#endif /* __mast_beam_optimization_h__ */
