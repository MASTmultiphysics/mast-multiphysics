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
#include "elasticity/piston_theory_boundary_condition.h"



MAST::PistonTheoryBoundaryCondition::
PistonTheoryBoundaryCondition(unsigned int order,
                              Real mach,
                              Real U_inf,
                              Real gamma,
                              Real rho,
                              const RealVectorX& vel_vec):
MAST::BoundaryConditionBase(MAST::PISTON_THEORY),
_order(order),
_mach(mach),
_U_inf(U_inf),
_gamma(gamma),
_rho(rho),
_vel_vec(vel_vec) {

    // scale the velocity vector to unit length
    _vel_vec.normalize();
}


MAST::PistonTheoryBoundaryCondition::~PistonTheoryBoundaryCondition() {
    
}



unsigned int
MAST::PistonTheoryBoundaryCondition::order() const {
    return _order;
}



Real
MAST::PistonTheoryBoundaryCondition::mach() const {
    return _mach;
}



Real
MAST::PistonTheoryBoundaryCondition::U_inf() const {
    return _U_inf;
}



void
MAST::PistonTheoryBoundaryCondition::set_U_inf(Real U_inf) {
    
    _U_inf = U_inf;
}



Real
MAST::PistonTheoryBoundaryCondition::rho() const {
    return _rho;
}



Real
MAST::PistonTheoryBoundaryCondition::gamma() const {
    return _gamma;
}



const RealVectorX&
MAST::PistonTheoryBoundaryCondition::vel_vec() const {
    return _vel_vec;
}


