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


// MAST includes
#include "fluid/pressure_function.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "numerics/utility.h"
#include "fluid/primitive_fluid_solution.h"
#include "fluid/small_disturbance_primitive_fluid_solution.h"
#include "fluid/flight_condition.h"
#include "base/nonlinear_system.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"
#include "libmesh/dof_map.h"


MAST::PressureFunction::
PressureFunction(MAST::SystemInitialization& sys,
                 MAST::FlightCondition&      flt):
MAST::FieldFunction<Real>("pressure"),
_if_cp            (false),
_ref_pressure     (0.),
_system           (sys),
_flt_cond         (flt) {
    
}




MAST::PressureFunction::~PressureFunction() {
    
}




void
MAST::PressureFunction::
init(const libMesh::NumericVector<Real>& steady_sol,
     const libMesh::NumericVector<Real>* small_dist_sol) {
    
    MAST::NonlinearSystem& sys = _system.system();
    
    // first initialize the solution to the given vector
    // steady state solution
    _sol.reset(libMesh::NumericVector<Real>::build(sys.comm()).release());
    _sol->init(steady_sol.size(), true, libMesh::SERIAL);
    
    // now localize the give solution to this objects's vector
    steady_sol.localize(*_sol);
    
    
    // if the mesh function has not been created so far, initialize it
    _sol_function.reset(new libMesh::MeshFunction(sys.get_equation_systems(),
                                                  *_sol,
                                                  sys.get_dof_map(),
                                                  _system.vars()));
    _sol_function->init();
    
    
    if (small_dist_sol) {
        
        // solution real part
        _dsol.reset(libMesh::NumericVector<Real>::build(sys.comm()).release());
        _dsol->init(steady_sol.size(), true, libMesh::SERIAL);
        
        small_dist_sol->localize(*_dsol);
        
        _dsol_function.reset(new libMesh::MeshFunction(sys.get_equation_systems(),
                                                       *_dsol,
                                                       sys.get_dof_map(),
                                                       _system.vars()));
        _dsol_function->init();
    }
    else {
        
        _dsol.reset();
        _dsol_function.reset();
    }
        
}





void
MAST::PressureFunction::
operator() (const libMesh::Point& p,
            const Real            time,
            Real                  &press) const {
    
    
    libmesh_assert(_sol_function.get()); // should be initialized before this call
    
    press  = 0.;
    
    
    // get the nonlinear and linearized solution
    DenseRealVector
    v;
    
    RealVectorX
    sol    = RealVectorX::Zero(_system.system().n_vars());
    
    
    // now the steady state function itself
    (*_sol_function)(p, 0., v);
    MAST::copy(sol, v);
    
    
    MAST::PrimitiveSolution                     p_sol;
    
    // now initialize the primitive variable contexts
    p_sol.init(dynamic_cast<MAST::ConservativeFluidSystemInitialization&>(_system).dim(),
               sol,
               _flt_cond.gas_property.cp,
               _flt_cond.gas_property.cv,
               false);
    
    if (_if_cp)
        press     = p_sol.c_pressure(_flt_cond.p0(), _flt_cond.q0());
    else {
        press     =  p_sol.p - _ref_pressure;
    }
}




void
MAST::PressureFunction::
perturbation(const libMesh::Point& p,
             const Real            t,
             Real                  &dpress) const {
    
    
    libmesh_assert(_sol_function.get()); // should be initialized before this call
    libmesh_assert(_dsol_function.get()); // should be initialized before this call
    
    dpress = 0.;
    
    
    // get the nonlinear and linearized solution
    DenseRealVector
    v;
    
    RealVectorX
    sol    = RealVectorX::Zero(_system.system().n_vars()),
    dsol   = RealVectorX::Zero(_system.system().n_vars());
    
    
    // first copy the real and imaginary solutions
    (*_dsol_function)(p, 0., v);
    MAST::copy(sol, v);
    dsol = sol;
    
    // now the steady state function itself
    (*_sol_function)(p, 0., v);
    MAST::copy(sol, v);
    
    
    MAST::PrimitiveSolution                     p_sol;
    SmallPerturbationPrimitiveSolution<Real> delta_p_sol;
    
    // now initialize the primitive variable contexts
    p_sol.init(dynamic_cast<MAST::ConservativeFluidSystemInitialization&>(_system).dim(),
               sol,
               _flt_cond.gas_property.cp,
               _flt_cond.gas_property.cv,
               false);
    delta_p_sol.init(p_sol,
                     dsol);
    
    if (_if_cp)
        dpress    = delta_p_sol.c_pressure(_flt_cond.q0());
    else
        dpress    =  delta_p_sol.dp;
}

