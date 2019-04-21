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

// C++ includes
#include <vector>

// MAST includes
#include "elasticity/structural_system.h"
#include "base/parameter.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/eigen_solver.h"
#include "libmesh/dof_map.h"


MAST::StructuralSystem::StructuralSystem(libMesh::EquationSystems& es,
                                         const std::string& name,
                                         const unsigned int number):
MAST::NonlinearSystem(es, name, number),
_iter(0),
_beta(0.),
_dl(0.),
_load_param(nullptr),
_min_p(0.),
_max_p(0.) {
    
}



MAST::StructuralSystem::~StructuralSystem() {
    
}



void
MAST::StructuralSystem::set_load_parameter(MAST::Parameter &param,
                                           Real min_p,
                                           Real max_p) {
    
    // make sure that it has not already been set
    libmesh_assert(!_load_param);
    
    _load_param = &param;
    _min_p      = min_p;
    _max_p      = max_p;
}



void
MAST::StructuralSystem::clear() {
    
    _load_param = nullptr;
    _beta       = 0.;
    _dl         = 0.;
    

    MAST::NonlinearSystem::clear();
}



void
MAST::StructuralSystem::solve(MAST::AssemblyElemOperations& elem_ops,
                              MAST::AssemblyBase& assembly) {

    // the system solves for both the structural displacement and
    // load update using a constraint definition
    //
    // the governing equations are
    //      r = 0    (structural residual)
    //      g = 0    (constraint equation)
    //
    //  the displacement and load steps are termed delta_x and delta_p. The
    //  Newton updates to these two quantities are termed dx and dp. Hence,
    //  each load increment seeks to solve for both delta_x and delta_p, while
    //  this nonlinear solution finds the incremental updates dx and dp so that
    //  the above two equations are satisfied.
    //
    //  for j as the load step and k as the N-R iterate,
    //  x_j = x_(j-1)  + delta_x_k
    //  p_j = p_(j-1)  + delta_p_k
    //
    //  delta_x_k = delta_x_(k-1) + dx_k
    //  delta_p_k = delta_p_(k-1) + dp_k
    //
    //  the linearized form is (where dx and dp are k_th N-R updates for j_th
    //  load step).
    //      dr/dx     dx +  dr/dp dp = -r
    //      (dg/dx)^T dx +  dg/dp dp = -g
    //
    //
    //   the system is solved as follows:
    //      dx1 = - inv(dr/dx) r
    //      dx2 = - inv(dr/dx) dr/dp
    //   so that
    //      dx  = dx1 + dp dx2, and
    // delta_x  = delta_x_(k-1) + dx
    //
    //   we use a quadratic constraint equation:
    //   g = delta_x^T delta_x + dl^2 = 0
    //
    //   we solve this in the following way:
    //     1. calculate dx1, dx2
    //     2. solve for dp as the solution of the quadratic constraint eq.
    //
    //   there is another way to solve this, by including the constraint in the
    //   N-R solution as follows. However, this will provide dp depending on
    //   the initial starting point for dp. Hence, we choose to analytically
    //   solve the quadratic equation.
    //
    //      (dg/dx)^T (dx1 + dp dx2)  + dg/dp dp  = -g
    //  or, ((dg/dx)^T dx2 + dg/dp) dp = -g - (dg/dx)^T dx1
    //  or, dp = - (g + (dg/dx)^T dx1) / ((dg/dx)^T dx2 + dg/dp)
    
    START_LOG("solve()", "StructuralSystem");

    //  the computations outlined above are repeated unless the
    //  nonlinear convergence criteria are satisfied.
    bool
    if_converged   = false;

    _dl = 1.;
    
    // these vectors will be used for solution updates here
    std::unique_ptr<libMesh::NumericVector<Real> >
    x_old        (this->solution->zero_clone().release()),  // to store solution after previous iteration
    delta_x0     (this->solution->zero_clone().release()),  // displacement update for this load step
    dx           (this->solution->zero_clone().release()),  // updates to delta_x
    dx_residual  (this->solution->zero_clone().release()),
    dx_load      (this->solution->zero_clone().release());

    // store the previous solution
    *x_old = *(this->solution);
    
    Real
    dp       = 0.,
    a0       = 0.,
    a1       = 0.,
    a2       = 0.,
    dp1      = 0.,
    dp2      = 0.,
    v1       = 0.;
    
    
    while (!if_converged) {
    
        // solve for dx due to residual. Set the maximum iterations to 1.

        // first solve the displacement update due to system residual
        // this solution will return x. We will use this to calculate
        // delta_x1 = x - x0
        //       dx = delta_x1 - delta_x0
        // where delta_x0 is the delta_x at the end of the previous iterate
        MAST::NonlinearSystem::solve(elem_ops, assembly);
        
        // if the total load steps have been taken to the max load, then
        // quit
        if (fabs((*_load_param)() - _max_p) < 1.0e-10) {
            if_converged = true;
            continue;
        }
        
        *dx = *(this->solution);
        dx->add(-1., *x_old);
        dx->add(-1., *delta_x0);
        
        // swap the solution for later use
        this->solution->swap(*dx_residual);
        
        
        // next solve the displacement update due to load update
        MAST::NonlinearSystem::sensitivity_solve(elem_ops, assembly, *_load_param);
        *dx_load = this->get_sensitivity_solution();
        
        // now, calculate the load update using the constraint definition.
        // We use a second order constraint function, which is written as a
        // quadratic equation
        // a0 dp^2 + a1 dp + a2 = 0
        a0   = dx_load->dot(*dx_load);
        a1   = 2.*(delta_x0->dot(*dx_load) + dx_residual->dot(*dx_load));
        v1   = -1.;  // some arbitrary value to trigger the while loop
        while (v1 < 0.) {
            
            a2   = (delta_x0->dot(*delta_x0)       +
                    2.*delta_x0->dot(*dx_residual) +
                    dx_residual->dot(*dx_residual) - _dl*_dl);
            v1   = a1*a1-4.*a0*a2;
            dp1  = 0.;
            dp2  = 0.;
            
            if (v1 >= 0.) {
                dp1  = (-a1 + sqrt(v1))/(2.*a0); // this is the first root
                dp2  = (-a1 - sqrt(v1))/(2.*a0); // this is the second root
            }
            else  {
                // reduce the value of dl and reevaluate the constraint
                _dl /= 2.;
            }
        }
        
        // now we need to select one of the two roots
        dp  = dp1;
        dp  = std::min(_max_p - (*_load_param)(), dp);
        
        // change the load  parameter
        (*_load_param)() += dp;
        
        // use the load update to evaluate the total displacement update
        *dx   =  *dx_residual;
        dx->add(dp, *dx_load);
        delta_x0->add(*dx);
        this->solution->add(1., *dx);
        
        // some preliminary convergence criteria
        if (dx->l2_norm() <= 1.e-6 ||
            fabs(dp)      <= 1.e-6)
            if_converged = true;
    }
    

    // Update the system after the solve
    this->update();

    STOP_LOG("solve()", "StructuralSystem");
}


