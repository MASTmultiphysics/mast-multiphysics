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


// MAST includes
#include "examples/fsi/plate_flag_flutter_solution/plate_flag_pressure_function.h"
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


MAST::PlateFlagPressureFunction::
PlateFlagPressureFunction(MAST::SystemInitialization& sys,
                         MAST::FlightCondition&      flt,
                         const Real th):
MAST::PressureFunction(sys, flt),
_flag_thickness        (th) {
    
}




MAST::PlateFlagPressureFunction::~PlateFlagPressureFunction() {
    
}





void
MAST::PlateFlagPressureFunction::
operator() (const libMesh::Point& p,
            const Real            time,
            Real                  &press) const {
    
    
    libmesh_assert(_sol_function.get()); // should be initialized before this call
    
    Real
    press_up  = 0.,
    press_lo  = 0.;
    press     = 0.;
    
    libMesh::Point
    pt;
    
    // get the nonlinear and linearized solution
    DenseRealVector
    v;
    
    RealVectorX
    sol    = RealVectorX::Zero(_system.system().n_vars());
    
    MAST::PrimitiveSolution                     p_sol;
    
    ///////////////////////////////////////////////////////////////
    //   pressure on the lower surface
    ///////////////////////////////////////////////////////////////
    
    pt = p;
    pt(2)  -=  0.5*_flag_thickness;
    
    // now the steady state function itself
    (*_sol_function)(pt, 0., v);
    MAST::copy(sol, v);
    
    // now initialize the primitive variable contexts
    p_sol.init(dynamic_cast<MAST::ConservativeFluidSystemInitialization&>(_system).dim(),
               sol,
               _flt_cond.gas_property.cp,
               _flt_cond.gas_property.cv,
               false);
    
    if (_if_cp)
        press_lo  = p_sol.c_pressure(_flt_cond.p0(), _flt_cond.q0());
    else
        press_lo  =  p_sol.p - _ref_pressure;
    
    ///////////////////////////////////////////////////////////////
    //   pressure on the upper surface
    ///////////////////////////////////////////////////////////////
    
    pt = p;
    pt(2)  +=  0.5*_flag_thickness;
    
    // now the steady state function itself
    (*_sol_function)(pt, 0., v);
    MAST::copy(sol, v);
    
    // now initialize the primitive variable contexts
    p_sol.zero();
    p_sol.init(dynamic_cast<MAST::ConservativeFluidSystemInitialization&>(_system).dim(),
               sol,
               _flt_cond.gas_property.cp,
               _flt_cond.gas_property.cv,
               false);
    
    if (_if_cp)
        press_up     = p_sol.c_pressure(_flt_cond.p0(), _flt_cond.q0());
    else
        press_up     =  p_sol.p - _ref_pressure;

    // the final pressure is the difference between the upper and lower
    // surfaces
    press = press_up - press_lo;
}




void
MAST::PlateFlagPressureFunction::
perturbation(const libMesh::Point& p,
             const Real            t,
             Real                  &dpress) const {
    
    
    libmesh_assert(_sol_function.get()); // should be initialized before this call
    libmesh_assert(_dsol_function.get()); // should be initialized before this call
    
    Real
    dpress_up  = 0.,
    dpress_lo  = 0.;
    dpress     = 0.;
    
    libMesh::Point
    pt;

    // get the nonlinear and linearized solution
    DenseRealVector
    v;
    
    RealVectorX
    sol    = RealVectorX::Zero(_system.system().n_vars()),
    dsol   = RealVectorX::Zero(_system.system().n_vars());
    
    MAST::PrimitiveSolution                     p_sol;
    SmallPerturbationPrimitiveSolution<Real> delta_p_sol;
    

    ///////////////////////////////////////////////////////////////
    //   pressure on the lower surface
    ///////////////////////////////////////////////////////////////
    
    pt = p;
    pt(2)  -=  0.5*_flag_thickness;

    // first copy the real and imaginary solutions
    (*_dsol_function)(pt, 0., v);
    MAST::copy(sol, v);
    dsol = sol;
    
    // now the steady state function itself
    (*_sol_function)(pt, 0., v);
    MAST::copy(sol, v);
    
    // now initialize the primitive variable contexts
    p_sol.init(dynamic_cast<MAST::ConservativeFluidSystemInitialization&>(_system).dim(),
               sol,
               _flt_cond.gas_property.cp,
               _flt_cond.gas_property.cv,
               false);
    delta_p_sol.init(p_sol,
                     dsol);
    
    if (_if_cp)
        dpress_lo    = delta_p_sol.c_pressure(_flt_cond.q0());
    else
        dpress_lo    =  delta_p_sol.dp;


    ///////////////////////////////////////////////////////////////
    //   pressure on the upper surface
    ///////////////////////////////////////////////////////////////
    
    pt = p;
    pt(2)  +=  0.5*_flag_thickness;
    
    // first copy the real and imaginary solutions
    (*_dsol_function)(pt, 0., v);
    MAST::copy(sol, v);
    dsol = sol;
    
    // now the steady state function itself
    (*_sol_function)(pt, 0., v);
    MAST::copy(sol, v);
    
    // now initialize the primitive variable contexts
    p_sol.zero();
    p_sol.init(dynamic_cast<MAST::ConservativeFluidSystemInitialization&>(_system).dim(),
               sol,
               _flt_cond.gas_property.cp,
               _flt_cond.gas_property.cv,
               false);
    delta_p_sol.zero();
    delta_p_sol.init(p_sol,
                     dsol);
    
    if (_if_cp)
        dpress_up    = delta_p_sol.c_pressure(_flt_cond.q0());
    else
        dpress_up    =  delta_p_sol.dp;
    
    // the final pressure is the difference between the upper and lower
    // surfaces
    dpress = dpress_up - dpress_lo;
}


