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

#ifndef __mast_fluid_elem_initialization_h__
#define __mast_fluid_elem_initialization_h__

// MAST includes
#include "base/nonlinear_system.h"
#include "base/transient_assembly.h"
#include "base/boundary_condition_base.h"
#include "fluid/conservative_fluid_element_base.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/conservative_fluid_transient_assembly.h"
#include "fluid/flight_condition.h"
#include "fluid/integrated_force_output.h"
#include "fluid/primitive_fluid_solution.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "mesh/geom_elem.h"
#include "base/test_comparisons.h"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_elem_type.h"    // ElemType
#include "libmesh/fe_type.h"           // FEFamily, Order
#include "libmesh/replicated_mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/mesh_generation.h"

extern libMesh::LibMeshInit* _libmesh_init;

struct BuildFluidElem {
    
    libMesh::LibMeshInit&                          _init;
    unsigned int                                   _dim;
    libMesh::UnstructuredMesh*                     _mesh;
    libMesh::EquationSystems*                      _eq_sys;
    MAST::NonlinearSystem*                         _sys;
    MAST::ConservativeFluidSystemInitialization*   _sys_init;
    MAST::ConservativeFluidDiscipline*             _discipline;
    MAST::FlightCondition*                         _flight_cond;
    MAST::TransientAssembly*                       _assembly;
    MAST::ConservativeFluidTransientAssemblyElemOperations* _elem_ops;
    MAST::FirstOrderNewmarkTransientSolver*        _transient_solver;
    MAST::GeomElem*                                _geom_elem;
    MAST::ConservativeFluidElementBase*            _fluid_elem;
    MAST::BoundaryConditionBase*                   _far_field_bc;
    MAST::BoundaryConditionBase*                   _slip_wall_bc;
    MAST::BoundaryConditionBase*                   _noslip_wall_bc;

    bool                                           _initialized;
    bool                                           _if_viscous;
    Real                                           _delta;
    libMesh::FEType                                _fetype;
    
    std::set<MAST::BoundaryConditionBase*>         _boundary_conditions;
    
    BuildFluidElem():
    _init           (*_libmesh_init),
    _dim            (0),
    _mesh           (nullptr),
    _eq_sys         (nullptr),
    _sys            (nullptr),
    _sys_init       (nullptr),
    _discipline     (nullptr),
    _flight_cond    (nullptr),
    _assembly       (nullptr),
    _elem_ops       (nullptr),
    _transient_solver(nullptr),
    _fluid_elem     (nullptr),
    _far_field_bc   (nullptr),
    _slip_wall_bc   (nullptr),
    _noslip_wall_bc (nullptr),
    _initialized    (false),
    _if_viscous     (false),
    _delta          (0.) {
        
    }
        
    void init(bool if_viscous) {
    
        // initialize the mesh. Details of parameters for each mesh are
        // described above.
        _mesh              = new libMesh::ReplicatedMesh(_init.comm());
        
        _dim               = 2;
        libMesh::MeshTools::Generation::build_square(*_mesh, 1, 1);
        
        // create equation system
        _eq_sys = new libMesh::EquationSystems(*_mesh);
        
        // add the system to be used for fluid analysis
        _sys = &(_eq_sys->add_system<MAST::NonlinearSystem>("fluid"));
        
        // create the discipline where boundary conditions will be stored
        _discipline = new MAST::ConservativeFluidDiscipline(*_eq_sys);
        
        // create system initialization object to add variables
        libMesh::Order    fe_order            = libMesh::FIRST;
        libMesh::FEFamily fe_family           = libMesh::LAGRANGE;
        
        _sys_init = new MAST::ConservativeFluidSystemInitialization(*_sys,
                                                                    _sys->name(),
                                                                    libMesh::FEType(fe_order, fe_family),
                                                                    _dim);
        
        _far_field_bc   = new MAST::BoundaryConditionBase(MAST::FAR_FIELD);
        _slip_wall_bc   = new MAST::BoundaryConditionBase(MAST::SLIP_WALL);
        _noslip_wall_bc = new MAST::BoundaryConditionBase(MAST::NO_SLIP_WALL);

        // set fluid properties
        _flight_cond = new MAST::FlightCondition;
        _flight_cond->flow_unit_vector(0)  = 1.;
        _flight_cond->flow_unit_vector(1)  = 1.;
        _flight_cond->flow_unit_vector(2)  = 0.;
        _flight_cond->mach                 = 0.5;
        _flight_cond->gas_property.cp      = 1003.;
        _flight_cond->gas_property.cv      = 716.;
        _flight_cond->gas_property.T       = 300.;
        _flight_cond->gas_property.rho     = 1.05;
        _flight_cond->gas_property.if_viscous = if_viscous;
        _flight_cond->init();
        
        // tell the discipline about the fluid values
        _discipline->set_flight_condition(*_flight_cond);
        
        // initialize the equation system
        _eq_sys->init();
        
        // not initialize the fluid element to be used for tests
        _assembly = new MAST::TransientAssembly;
        _elem_ops = new MAST::ConservativeFluidTransientAssemblyElemOperations;
        _transient_solver = new MAST::FirstOrderNewmarkTransientSolver;
        
        _assembly->set_discipline_and_system(*_discipline, *_sys_init);
        _elem_ops->set_discipline_and_system(*_discipline, *_sys_init);
        _transient_solver->set_discipline_and_system(*_discipline, *_sys_init);
        _transient_solver->set_elem_operation_object(*_elem_ops);
        _geom_elem = new MAST::GeomElem;
        _geom_elem->init(**_mesh->elements_begin(), *_sys_init);
        _fluid_elem = new MAST::ConservativeFluidElementBase(*_sys_init,
                                                             *_geom_elem,
                                                             *_flight_cond);
        
        _initialized = true;
    }
    
    
    ~BuildFluidElem() {
        
        if (_initialized) {
            
            delete _assembly;
            delete _elem_ops;
            delete _transient_solver;
            delete _fluid_elem;
            delete _geom_elem;
            
            delete _far_field_bc;
            delete _slip_wall_bc;
            delete _noslip_wall_bc;
            
            delete _eq_sys;
            delete _mesh;
            
            delete _discipline;
            delete _sys_init;
            delete _flight_cond;
        }
    }
    
    
    void init_solution_for_elem(RealVectorX& s) {
        
        libmesh_assert(_initialized);
        
        unsigned int
        n_shape = _sys->n_dofs()/(_dim+2);
        
        // random vector in [-1, 1] used to perturb the state vector
        // to invoke gradients, leading to stresses and flux terms
        RealVectorX
        rand = RealVectorX::Random(n_shape);
        
        s = RealVectorX::Zero(_sys->n_dofs());
        
        for (unsigned int i=0; i<n_shape; i++) {
            
            s(i)                  = _flight_cond->rho()    * (1. + _delta * rand(i));
            s(n_shape+i)          = _flight_cond->rho_u1() * (1. + _delta * rand(i));
            s(2*n_shape+i)        = _flight_cond->rho_u2() * (1. + _delta * rand(i));
            if (_dim > 2)
                s(3*n_shape+i)    = _flight_cond->rho_u3() * (1. + _delta * rand(i));
            s((_dim+1)*n_shape+i) = _flight_cond->rho_e()  * (1. + _delta * rand(i));
        }
    }



    void init_primitive_sol(MAST::PrimitiveSolution& s) {
     
        libmesh_assert(_initialized);
        
        RealVectorX
        c_sol = RealVectorX::Zero(_dim+2);
        
        c_sol(0)                  = _flight_cond->rho();
        c_sol(1)                  = _flight_cond->rho_u1();
        c_sol(2)                  = _flight_cond->rho_u2();
        if (_dim > 2)
            c_sol(3)    = _flight_cond->rho_u3();
        c_sol(_dim+1) = _flight_cond->rho_e();
        
        s.zero();
        s.init(_dim, c_sol,
               _flight_cond->gas_property.cp,
               _flight_cond->gas_property.cv,
               _flight_cond->gas_property.if_viscous);
    }
    
    template <typename ValType>
    bool check_jacobian(ValType& val) {
        
        libmesh_assert(_initialized);

        Real
        delta = 0.;
        
        unsigned int
        n = _sys->n_dofs();
        
        RealMatrixX
        jac0    = RealMatrixX::Zero(n, n),
        jac_fd  = RealMatrixX::Zero(n, n);
        
        RealVectorX
        v_low   = RealVectorX::Zero(n),
        v_hi    = RealVectorX::Zero(n),
        v0      = RealVectorX::Zero(n),
        x_low   = RealVectorX::Zero(n),
        x_hi    = RealVectorX::Zero(n),
        x0      = RealVectorX::Zero(n),
        f0      = RealVectorX::Zero(n),
        f_low   = RealVectorX::Zero(n),
        f_hi    = RealVectorX::Zero(n);
        
        init_solution_for_elem(x0);
        
        // velocity is set to assuming a linear variation of the state from zero
        // over dt = 1.e-2
        v0 = x0/1.e-2;
        
        _fluid_elem->set_solution(x0);
        _fluid_elem->set_velocity(v0);
        
        val.compute(true, f0, jac0);
        
        for (unsigned int i=0; i<n; i++) {
            
            x_low = x0;
            x_hi  = x0;
            v_low = v0;
            v_hi  = v0;
            
            if (!val.jac_xdot) {
                if (fabs(x0(i)) > 0.)
                    delta = x0(i)*val.frac;
                else
                    delta = val.delta;
                
                x_low(i) -= delta;
                x_hi(i)  += delta;
            }
            else {
                if (fabs(v0(i)) > 0.)
                    delta = v0(i)*val.frac;
                else
                    delta = val.delta;
                
                v_low(i) -= delta;
                v_hi(i)  += delta;
            }
            
            
            // get the new residual
            _fluid_elem->set_solution(x_low);
            _fluid_elem->set_velocity(v_low);
            f_low.setZero();
            val.compute(false, f_low, jac0);
            
            
            // get the new residual
            _fluid_elem->set_solution(x_hi);
            _fluid_elem->set_velocity(v_hi);
            f_hi.setZero();
            val.compute(false, f_hi, jac0);
            
            // set the i^th column of the finite-differenced Jacobian
            jac_fd.col(i) = (f_hi-f_low)/2./delta;
        }
        
        return MAST::compare_matrix(jac_fd, jac0, val.tol);
    }
    
};

#endif // __mast_fluid_elem_initialization_h__
