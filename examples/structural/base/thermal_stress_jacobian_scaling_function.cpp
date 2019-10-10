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
#include "examples/structural/base/thermal_stress_jacobian_scaling_function.h"
#include "base/nonlinear_implicit_assembly.h"


MAST::Examples::ThermalJacobianScaling::ThermalJacobianScaling():
MAST::FieldFunction<Real>("thermal_jacobian_scaling"),
_assembly     (nullptr),
_enable       (false),
_accel_factor (2.)    {
    
}


void
MAST::Examples::ThermalJacobianScaling::set_assembly(MAST::NonlinearImplicitAssembly& assembly) {

    libmesh_assert(!_assembly);
    _assembly = &assembly;
}


void
MAST::Examples::ThermalJacobianScaling::clear_assembly() {

    _assembly = nullptr;
}


void
MAST::Examples::ThermalJacobianScaling::set_enable(bool f) {
    
    _enable = f;
}


void
MAST::Examples::ThermalJacobianScaling::set_acceleration_factor(Real f) {

    libmesh_assert_greater(f, 0.);
    
    _accel_factor = f;
}



void
MAST::Examples::ThermalJacobianScaling::operator() (Real& val) const {
    
    libmesh_assert(_assembly);
    
    if (!_enable) {
        
        val  =  1.;
        return;
    }
    
    Real
    low      = 1.e-2,
    res      = _assembly->res_l2_norm(),
    max_res  = _assembly->first_iter_res_l2_norm(),
    max_log  = 0.,
    log      = res>0.?std::log10(res):0.;
    
    // scaling is at low value for highest norm, otherwise
    // it gets scaled
    if (max_res < 0.) val      =  low;
    else {
        
        if (max_res > 0.)
            max_log = std::log10(max_res);
        val      = std::fmax(std::pow(1-res/max_res, _accel_factor), low);
        //libMesh::out << log << " " << max_log << " " << val << std::endl;
    }
}

