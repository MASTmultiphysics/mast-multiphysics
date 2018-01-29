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
#include "fluid/surface_integrated_pressure_output.h"



MAST::SurfaceIntegratedPressureOutput::
SurfaceIntegratedPressureOutput(MAST::SurfaceIntegratedPressureOutput::OutputMode o,
                                const RealVectorX& n_vec):
MAST::OutputAssemblyElemOperations(),
_mode(o),
_n_vec(n_vec) {

    // scale the vector if needed
    if (o == MAST::SurfaceIntegratedPressureOutput::OutputMode::UNIT_VEC &&
        n_vec.norm() >= 0.)
        _n_vec /= _n_vec.norm();
}


MAST::SurfaceIntegratedPressureOutput::
~SurfaceIntegratedPressureOutput() {
    
}


void
MAST::SurfaceIntegratedPressureOutput::clear() {
    
    _load.setZero();
    _load_sensitivity.clear();
    _dload_dX.setZero();
}



void
MAST::SurfaceIntegratedPressureOutput::
set_output_mode(MAST::SurfaceIntegratedPressureOutput::OutputMode o,
                const RealVectorX* n_vec) {

    _mode = o;
    
    // scale the vector if needed
    if (o == MAST::SurfaceIntegratedPressureOutput::OutputMode::UNIT_VEC) {
        
        libmesh_assert(n_vec);
        _n_vec = *n_vec;
        _n_vec /= _n_vec.norm();
    }
}

    

Real
MAST::SurfaceIntegratedPressureOutput::
value() const {
 
    switch (_mode) {
        case MAST::SurfaceIntegratedPressureOutput::L2_NORM: {
            
            return _load.norm();
        }
            break;

        case MAST::SurfaceIntegratedPressureOutput::UNIT_VEC: {
            
            return _load.dot(_n_vec);
        }
            break;

        default:
            libmesh_assert(false); // should not get here
            break;
    }
}
        
        

Real
MAST::SurfaceIntegratedPressureOutput::
sensitivity(const MAST::FunctionBase& f) const {
    
    // get the sensitivity of the load wrt the specified parameter
    std::map<const MAST::FunctionBase*, RealVectorX>::const_iterator
    it = _load_sensitivity.find(&f);
    
    libmesh_assert(it != _load_sensitivity.end());
    
    switch (_mode) {
        case MAST::SurfaceIntegratedPressureOutput::L2_NORM: {
            
            return (_load.dot(it->second)) / _load.norm();
        }
            break;
            
        case MAST::SurfaceIntegratedPressureOutput::UNIT_VEC: {
            
            return it->second.dot(_n_vec);
        }
            break;
            
        default:
            libmesh_assert(false); // should not get here
            break;
    }
}
        
        

RealVectorX
MAST::SurfaceIntegratedPressureOutput::derivative() const {

    switch (_mode) {
        case MAST::SurfaceIntegratedPressureOutput::L2_NORM: {
            
            return (_load.transpose()*_dload_dX) / _load.norm();
        }
            break;
            
        case MAST::SurfaceIntegratedPressureOutput::UNIT_VEC: {
            
            return _n_vec.transpose() * _dload_dX;
        }
            break;
            
        default:
            libmesh_assert(false); // should not get here
            break;
    }
}
        
