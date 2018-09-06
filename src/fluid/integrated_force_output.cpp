/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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
#include "fluid/integrated_force_output.h"
#include "fluid/conservative_fluid_element_base.h"
#include "fluid/conservative_fluid_discipline.h"
#include "base/assembly_base.h"
#include "base/elem_base.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"


// libMesh includes
#include "libmesh/boundary_info.h"


MAST::IntegratedForceOutput::IntegratedForceOutput(const RealVectorX& nvec):
MAST::OutputAssemblyElemOperations(),
_n_vec   (nvec),
_force   (0.) {
    
    _n_vec /= _n_vec.norm();
}



MAST::IntegratedForceOutput::~IntegratedForceOutput()  {
    
}



void
MAST::IntegratedForceOutput::init(const libMesh::Elem& elem) {
    
    libmesh_assert(!_physics_elem);
    libmesh_assert(_system);
    libmesh_assert(_assembly);
    
    const MAST::FlightCondition& p =
    dynamic_cast<MAST::ConservativeFluidDiscipline&>
    (_assembly->discipline()).flight_condition();
    
    _physics_elem =
    new MAST::ConservativeFluidElementBase(*_system, *_assembly, elem, p);
}



void
MAST::IntegratedForceOutput::zero_for_analysis()  {
    
    _force = 0.;
}


void
MAST::IntegratedForceOutput::zero_for_sensitivity() {
    
    // nothing to be done here
}


Real
MAST::IntegratedForceOutput::output_total()  {

    Real val = _force;
    _system->system().comm().sum(val);
    return val;
}



void
MAST::IntegratedForceOutput::evaluate() {

    libmesh_assert(_physics_elem);

    MAST::ConservativeFluidElementBase& e =
    dynamic_cast<MAST::ConservativeFluidElementBase&>(*_physics_elem);

    const libMesh::Elem&
    elem = _physics_elem->elem();
    
    RealVectorX
    f  = RealVectorX::Zero(3);
    
    for (unsigned short int n=0; n<elem.n_sides(); n++)
        if (this->if_evaluate_for_boundary(elem, n)) {
            
            e.side_integrated_force(n, f);
            _force += f.dot(_n_vec);
        }
}


void
MAST::IntegratedForceOutput::output_derivative_for_elem(RealVectorX& dq_dX) {
    
    libmesh_assert(_physics_elem);
    
    MAST::ConservativeFluidElementBase& e =
    dynamic_cast<MAST::ConservativeFluidElementBase&>(*_physics_elem);
    
    const libMesh::Elem&
    elem = _physics_elem->elem();
    
    RealVectorX
    f    = RealVectorX::Zero(3);
    
    RealMatrixX
    dfdX = RealMatrixX::Zero(3, dq_dX.size());
    
    for (unsigned short int n=0; n<elem.n_sides(); n++)
        if (this->if_evaluate_for_boundary(elem, n)) {
            
            e.side_integrated_force(n, f, &dfdX);
            dq_dX += _n_vec.transpose() * dfdX;
        }
}



