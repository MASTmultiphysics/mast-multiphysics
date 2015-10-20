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

// C++ includes
#include <iomanip>


// MAST includes
#include "fluid/small_disturbance_primitive_fluid_solution.h"
#include "fluid/primitive_fluid_solution.h"


template <typename ValType>
MAST::SmallPerturbationPrimitiveSolution<ValType>::
SmallPerturbationPrimitiveSolution() {
    
    this->zero();
}


template <typename ValType>
void MAST::SmallPerturbationPrimitiveSolution<ValType>::zero() {
    
    perturb_primitive_sol.setZero();
    
    drho          = 0.;
    du1           = 0.;
    du2           = 0.;
    du3           = 0.;
    dT            = 0.;
    dp            = 0.;
    da            = 0.;
    de_tot        = 0.;
    dk            = 0.;
    dentropy      = 0.;
    dmach         = 0.;
    primitive_sol = NULL;
}


template <typename ValType>
void
MAST::SmallPerturbationPrimitiveSolution<ValType>::
init(const MAST::PrimitiveSolution& sol,
     const typename VectorType<ValType>::return_type& delta_sol) {
    
    primitive_sol = &sol;
    
    const unsigned int
    n1        = sol.dimension+2;
    
    const
    Real
    R         = sol.cp-sol.cv,
    gamma     = sol.cp/sol.cv;
    
    perturb_primitive_sol.resize(n1);
    
    drho = delta_sol(0);
    perturb_primitive_sol(0) = drho;
    
    du1 = (delta_sol(1) - drho * sol.u1)/sol.rho;
    
    perturb_primitive_sol(1) = du1;
    
    dk = sol.u1*du1;
    
    if (sol.dimension > 1)
    {
        du2 = (delta_sol(2) - drho * sol.u2)/sol.rho;
        perturb_primitive_sol(2) = du2;
        dk += sol.u2*du2;
    }
    
    if (sol.dimension > 2)
    {
        du3 = (delta_sol(3) - drho * sol.u3)/sol.rho;
        perturb_primitive_sol(3) = du3;
        dk += sol.u3*du3;
    }
    
    de_tot = (delta_sol(n1-1) - drho * sol.e_tot)/sol.rho;
    
    dT = (de_tot - dk)/sol.cv;
    perturb_primitive_sol(n1-1) = dT;
    
    dp = R*(dT*sol.rho + sol.T*drho);
    da = 0.5*sqrt(gamma*R/sol.T)*dT;
    dmach =  dk/sqrt(2.0*sol.k)/sol.a - sqrt(2.*sol.k)/pow(sol.a,2) * da;
    dentropy = (dp/pow(sol.rho,gamma) - gamma*sol.p/pow(sol.rho,gamma+1.)*drho)
    / (sol.p/pow(sol.rho,gamma)) ;
}


template <typename ValType>
void
MAST::SmallPerturbationPrimitiveSolution<ValType>::
print(std::ostream& out) const {
    
    out
    << "Small Perturbation Primitive Solution:" << std::endl
    << perturb_primitive_sol << std::endl
    << std::setw(15) <<  " drho: " << drho << std::endl
    << std::setw(15) <<  " du1: " << du1 << std::endl
    << std::setw(15) <<  " du2: " << du2 << std::endl
    << std::setw(15) <<  " du3: " << du3 << std::endl
    << std::setw(15) <<  " dmach: " << dmach << std::endl
    << std::setw(15) <<  " da: " << da << std::endl
    << std::setw(15) <<  " dT: " << dT << std::endl
    << std::setw(15) <<  " dp: " << dp << std::endl
    << std::setw(15) <<  " de_tot: " << de_tot << std::endl
    << std::setw(15) <<  " dk: " << dk << std::endl
    << std::setw(15) <<  " dentropy: " << dentropy << std::endl << std::endl;
}


template <typename ValType>
ValType
MAST::SmallPerturbationPrimitiveSolution<ValType>::
c_pressure(const Real q0) const {
    
    return dp/q0;
}


template <typename ValType>
void
MAST::SmallPerturbationPrimitiveSolution<ValType>::
get_duvec(typename VectorType<ValType>::return_type& du) const {
    
    du.setZero();
    
    switch (primitive_sol->dimension) {
        case 3:
            du(2) = du3;
        case 2:
            du(1) = du2;
        case 1:
            du(0) = du1;
    }
}



// explicit instantiations for real and complex type
template class MAST::SmallPerturbationPrimitiveSolution<Real>;
template class MAST::SmallPerturbationPrimitiveSolution<Complex>;

