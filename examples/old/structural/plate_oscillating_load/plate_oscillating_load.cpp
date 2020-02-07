/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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
#include "examples/structural/plate_oscillating_load/plate_oscillating_load.h"
#include "examples/base/input_wrapper.h"
#include "base/field_function_base.h"
#include "base/parameter.h"
#include "base/boundary_condition_base.h"
#include "base/physics_discipline_base.h"

namespace MAST {
    namespace Examples {
        
        // this is a function definition for the oscillating pressure
        class OscillatingDistributedLoad:
        public MAST::FieldFunction<Real> {
            
        public:
            
            /*!
             *   \p p is the distributed load and \p f is the circular frequency
             */
            OscillatingDistributedLoad(MAST::Parameter& p,
                                       MAST::Parameter& f):
            MAST::FieldFunction<Real>("pressure"),
            _p(p),
            _f(f) {
                
                _functions.insert(&p);
                _functions.insert(&f);
            }
            
            virtual ~OscillatingDistributedLoad() { }
            
            
            /*!
             *    calculates the value of the function at the specified point,
             *    \p p, and time, \p t, and returns it in \p v.
             */
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& v) const {
                
                v = _p() * sin(_f()*t);
            }
            
            
            /*!
             *    calculates the value of the function at the specified point,
             *    \p p, and time, \p t, and returns it in \p v.
             */
            virtual void derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& v) const {
                
                Real
                dp = 0.,
                df = 0.;
                
                // if the sensitivity parameter is the load parameter itself,
                // then the sensitivity parameter will be nonzero.
                if (_p.depends_on(f)) dp = 1.;
                if (_f.depends_on(f)) df = 1.;
                
                v = dp * sin(_f()*t) + _p() * _f() * df * cos(_f() * t);
            }
            
        protected:
            
            MAST::Parameter& _p;
            
            MAST::Parameter& _f;
            
        };
    }
}


MAST::Examples::PlateOscillatingLoad::
PlateOscillatingLoad(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::StructuralExample2D(comm_in) {
    
}



void
MAST::Examples::PlateOscillatingLoad::_init_loads() {
    
    Real
    p_val    =  (*_input)(_prefix+"pressure",                     "pressure on plate surface",     2.e4),
    freq_val =  (*_input)(   _prefix+"omega",  "frequency of oscillation of pressure (rad/s)",     1.e2);
    
    MAST::Parameter
    *press   = new MAST::Parameter( "p",         p_val),
    *freq    = new MAST::Parameter( "omega",  freq_val);
    
    MAST::Examples::OscillatingDistributedLoad
    *press_f         = new MAST::Examples::OscillatingDistributedLoad(*press, *freq);
    
    // initialize the load
    MAST::BoundaryConditionBase
    *p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
    
    p_load->add(*press_f);
    _discipline->add_volume_load(0, *p_load);
    
    
    this->add_parameter(*press);
    this->add_parameter(*freq);
    this->register_field_function(*press_f);
    this->register_loading(*p_load);
}

