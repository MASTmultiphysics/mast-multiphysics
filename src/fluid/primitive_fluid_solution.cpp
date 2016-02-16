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


// C++ includes
#include <iomanip>

// MAST includes
#include "fluid/primitive_fluid_solution.h"


MAST::PrimitiveSolution::PrimitiveSolution()
{
    this->zero();
}


void
MAST::PrimitiveSolution::zero()
{
    this->primitive_sol.setZero();
    dimension = 0;
    cp = 0.;
    cv = 0.;
    rho = 0.;
    u1 = 0.;
    u2 = 0.;
    u3 = 0.;
    T = 0.;
    p = 0.;
    a = 0.;
    e_tot = 0.;
    k = 0.;
    entropy = 0.;
    mach = 0.;
    // viscous values
    Pr = 0.72;
    k_thermal = 0.;
    mu = 0.;
    lambda = 0.;
}



void
MAST::PrimitiveSolution::init(const unsigned int dim,
                              const RealVectorX &conservative_sol,
                              const Real cp_val,
                              const Real cv_val,
                              bool if_viscous)
{
    dimension = dim;
    const unsigned int n1 = dim+2;
    cp = cp_val; cv = cv_val;
    const Real R = cp-cv, gamma = cp/cv;
    primitive_sol.resize(n1);
    
    rho = conservative_sol(0);
    primitive_sol(0) = rho;
    
    u1 = conservative_sol(1)/rho;
    primitive_sol(1) = u1;
    k = u1*u1;
    
    if (dim > 1)
    {
        u2 = conservative_sol(2)/rho;
        primitive_sol(2) = u2;
        k += u2*u2;
    }
    
    if (dim > 2)
    {
        u3 = conservative_sol(3)/rho;
        primitive_sol(3) = u3;
        k += u3*u3;
    }
    k *= 0.5;
    
    e_tot = conservative_sol(n1-1)/rho; // cv*T+k;
    
    T = (e_tot - k)/cv;
    primitive_sol(n1-1) = T;
    
    p = R*T*rho;
    a = sqrt(gamma*R*T);
    mach = sqrt(2.0*k)/a;
    entropy = log(p/pow(rho,gamma));
    
    // viscous quantities
    if (if_viscous)
    {
        mu = 1.458e-6 * pow(T, 1.5)/(T+110.4);
        lambda = -2./3.*mu;
        k_thermal = mu*cp/Pr;
    }
}



Real
MAST::PrimitiveSolution::c_pressure(const Real p0,
                                    const Real q0) const
{
    return (p-p0)/q0;
}



void
MAST::PrimitiveSolution::get_uvec(RealVectorX &u) const
{
    u.setZero();
    
    switch (dimension) {
        case 3:
            u(2) = u3;
        case 2:
            u(1) = u2;
        case 1:
            u(0) = u1;
    }
}



void
MAST::PrimitiveSolution::print(std::ostream& out) const
{
    out
    << "Primitive Solution:" << std::endl
    << primitive_sol << std::endl
    << std::setw(15) <<  " rho: " << rho << std::endl
    << std::setw(15) <<  " u1: " << u1 << std::endl
    << std::setw(15) <<  " u2: " << u2 << std::endl
    << std::setw(15) <<  " u3: " << u3 << std::endl
    << std::setw(15) <<  " mach: " << mach << std::endl
    << std::setw(15) <<  " a: " << a << std::endl
    << std::setw(15) <<  " T: " << T << std::endl
    << std::setw(15) <<  " p: " << p << std::endl
    << std::setw(15) <<  " e_tot: " << e_tot << std::endl
    << std::setw(15) <<  " k: " << k << std::endl
    << std::setw(15) <<  " entropy: " << entropy << std::endl
    << std::setw(15) <<  " Pr: " << Pr << std::endl
    << std::setw(15) <<  " mu: " << mu << std::endl
    << std::setw(15) <<  " lambda: " << lambda << std::endl
    << std::setw(15) <<  " k_thermal: " << k_thermal << std::endl << std::endl;
    
}


