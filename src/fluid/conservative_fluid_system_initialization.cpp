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
#include "fluid/conservative_fluid_system_initialization.h"


MAST::ConservativeFluidSystemInitialization::
ConservativeFluidSystemInitialization(libMesh::System& sys,
                                      const std::string& prefix,
                                      const libMesh::FEType& fe_type,
                                      const unsigned int dim):
MAST::SystemInitialization(sys, prefix) {
    
    _vars.resize(dim+2);
    
    std::string nm = prefix + "_rho";
    _vars[0] = sys.add_variable(nm, fe_type);
    
    nm = prefix + "_rhoux";
    _vars[1] = sys.add_variable(nm, fe_type);
    
    if (dim > 1) {
        nm = prefix + "_rhouy";
        _vars[2] = sys.add_variable(nm, fe_type);
    }

    if (dim > 2) {
        nm = prefix + "_rhouz";
        _vars[3] = sys.add_variable(nm, fe_type);
    }

    
    nm = prefix + "_rhoe";
    _vars[dim+2] = sys.add_variable(nm, fe_type);
}


MAST::ConservativeFluidSystemInitialization::
~ConservativeFluidSystemInitialization() {
    
}
