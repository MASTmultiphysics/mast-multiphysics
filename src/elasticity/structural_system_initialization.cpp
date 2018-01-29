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
#include "elasticity/structural_system_initialization.h"
#include "base/nonlinear_system.h"


// libMesh includes
#include "libmesh/equation_systems.h"


MAST::StructuralSystemInitialization::
StructuralSystemInitialization(MAST::NonlinearSystem& sys,
                               const std::string& prefix,
                               const libMesh::FEType& fe_type):
MAST::SystemInitialization(sys, prefix),
_stress_output_sys(nullptr) {
    
    _vars.resize(6);
    
    std::string nm = prefix + "_ux";
    _vars[0] = sys.add_variable(nm, fe_type);

    nm = prefix + "_uy";
    _vars[1] = sys.add_variable(nm, fe_type);

    nm = prefix + "_uz";
    _vars[2] = sys.add_variable(nm, fe_type);

    nm = prefix + "_tx";
    _vars[3] = sys.add_variable(nm, fe_type);

    nm = prefix + "_ty";
    _vars[4] = sys.add_variable(nm, fe_type);

    nm = prefix + "_tz";
    _vars[5] = sys.add_variable(nm, fe_type);

    
    // now initialize the stress system for output of stress data
    _stress_output_sys  =
    &(sys.get_equation_systems().add_system<libMesh::ExplicitSystem>("StressOutput"));
    
    libMesh::FEType
    fetype(libMesh::CONSTANT, libMesh::MONOMIAL); // constant value per element
    
    _stress_vars.resize(13);
    
    _stress_vars[0]  = _stress_output_sys->add_variable("epsilon-xx", fetype);
    _stress_vars[1]  = _stress_output_sys->add_variable("epsilon-yy", fetype);
    _stress_vars[2]  = _stress_output_sys->add_variable("epsilon-zz", fetype);
    _stress_vars[3]  = _stress_output_sys->add_variable("epsilon-xy", fetype);
    _stress_vars[4]  = _stress_output_sys->add_variable("epsilon-yz", fetype);
    _stress_vars[5]  = _stress_output_sys->add_variable("epsilon-zx", fetype);
    
    _stress_vars[6]  = _stress_output_sys->add_variable("sigma-xx",    fetype);
    _stress_vars[7]  = _stress_output_sys->add_variable("sigma-yy",    fetype);
    _stress_vars[8]  = _stress_output_sys->add_variable("sigma-zz",    fetype);
    _stress_vars[9]  = _stress_output_sys->add_variable("sigma-xy",    fetype);
    _stress_vars[10] = _stress_output_sys->add_variable("sigma-yz",    fetype);
    _stress_vars[11] = _stress_output_sys->add_variable("sigma-zx",    fetype);
    
    _stress_vars[12] = _stress_output_sys->add_variable("sigma-vm",    fetype);
    
}


MAST::StructuralSystemInitialization::~StructuralSystemInitialization() {
    
}
