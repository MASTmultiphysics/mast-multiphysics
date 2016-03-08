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


// MAST includes
#include "fluid/small_disturbance_pressure_function.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "numerics/utility.h"
#include "fluid/primitive_fluid_solution.h"
#include "fluid/small_disturbance_primitive_fluid_solution.h"
#include "fluid/flight_condition.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"
#include "libmesh/dof_map.h"


MAST::SmallDisturbancePressureFunction::
SmallDisturbancePressureFunction(MAST::SystemInitialization& sys,
                                 MAST::FlightCondition&      flt):
MAST::SmallDisturbancePressure(),
_system(sys),
_flt_cond(flt) {
    
}




MAST::SmallDisturbancePressureFunction::~SmallDisturbancePressureFunction() {
    
}




void
MAST::SmallDisturbancePressureFunction::
init(libMesh::NumericVector<Real>& steady_sol,
     libMesh::NumericVector<Real>& small_dist_sol_real,
     libMesh::NumericVector<Real>& small_dist_sol_imag) {
    
    libMesh::System& sys = _system.system();
    
    // first initialize the solution to the given vector
    // steady state solution
    _sol.reset(libMesh::NumericVector<Real>::build(sys.comm()).release());
    _sol->init(sys.n_dofs(),
               sys.n_local_dofs(),
               sys.get_dof_map().get_send_list(),
               false,
               libMesh::GHOSTED);
    
    // solution real part
    _dsol_real.reset(libMesh::NumericVector<Real>::build(sys.comm()).release());
    _dsol_real->init(sys.n_dofs(),
                     sys.n_local_dofs(),
                     sys.get_dof_map().get_send_list(),
                     false,
                     libMesh::GHOSTED);
    
    // solution complex part
    _dsol_imag.reset(libMesh::NumericVector<Real>::build(sys.comm()).release());
    _dsol_imag->init(sys.n_dofs(),
                     sys.n_local_dofs(),
                     sys.get_dof_map().get_send_list(),
                     false,
                     libMesh::GHOSTED);
    
    // now localize the give solution to this objects's vector
    steady_sol.localize(*_sol);
    small_dist_sol_real.localize(*_dsol_real);
    small_dist_sol_imag.localize(*_dsol_imag);
    
    
    // if the mesh function has not been created so far, initialize it
    _function.reset(new libMesh::MeshFunction(_system.system().get_equation_systems(),
                                              *_sol,
                                              _system.system().get_dof_map(),
                                              _system.vars()));
    _function->init();
}





void
MAST::SmallDisturbancePressureFunction::
freq_domain_pressure(const libMesh::Point& p,
                     const bool if_cp,
                     Real&    press,
                     Complex& dpress) {
    
    
    libmesh_assert(_function.get()); // should be initialized before this call
    
    press  = 0.;
    dpress = 0.;
    
    
    // get the nonlinear and linearized solution
    DenseRealVector
    v;
    
    RealVectorX
    sol    = RealVectorX::Zero(_system.system().n_vars());
    ComplexVectorX
    dsol   = ComplexVectorX::Zero(_system.system().n_vars());
    
    
    // first copy the real and imaginary solutions
    _sol->swap(*_dsol_real);
    (*_function)(p, 0., v);
    _sol->swap(*_dsol_real);
    MAST::copy(sol, v);
    dsol.real() = sol;
    
    
    // now the imaginary part
    _sol->swap(*_dsol_imag);
    (*_function)(p, 0., v);
    _sol->swap(*_dsol_imag);
    MAST::copy(sol, v);
    dsol.imag() = sol;
    
    
    // now the steady state function itself
    (*_function)(p, 0., v);
    MAST::copy(sol, v);
    
    
    MAST::PrimitiveSolution                     p_sol;
    SmallPerturbationPrimitiveSolution<Complex> delta_p_sol;
    
    // now initialize the primitive variable contexts
    p_sol.init(dynamic_cast<MAST::ConservativeFluidSystemInitialization&>(_system).dim(),
               sol,
               _flt_cond.gas_property.cp,
               _flt_cond.gas_property.cv,
               false);
    delta_p_sol.init(p_sol,
                     dsol);
    
    if (if_cp) {
        
        press     = p_sol.c_pressure(_flt_cond.p0(),
                                  _flt_cond.q0());
        dpress    = delta_p_sol.c_pressure(_flt_cond.q0());
    }
    else {
        
        press     =  p_sol.p;
        dpress    =  delta_p_sol.dp;
    }
}





