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

#ifndef __mast_fluid_condition_h__
#define __mast_fluid_condition_h__

// MAST includes
#include "base/mast_data_types.h"
#include "fluid/gas_property.h"
#include "base/field_function_base.h"


// Eigen includes
#include "Eigen/Core"

namespace MAST {
    
    class FlightCondition {
        
    public:
        FlightCondition():
        flow_unit_vector       (RealVectorX::Zero(3)),
        ref_chord              (0.)
        enable_shock_capturing (true)
        {}
        
        virtual ~FlightCondition()
        {}

        /*!
         *  flag to turn on/off artificial dissipation for
         *  shock capturing terms in fluid flow analysis. This is
         *  \p true by default.
         */
        bool enable_shock_capturing;
        
        /*!
         *   direction along which flow is defined
         */
        RealVectorX  flow_unit_vector;
        
        /*!
         *   Flight Mach number
         */
        Real mach;
        
        /*!
         *  Velocity magnitude, whose direction is evaluated from the Euler angles.
         */
        Real velocity_magnitude() const {return mach * gas_property.a;}
        
        /*!
         *   Ambient air properties
         */
        GasProperty gas_property;
        
        /*!
         *   reference chord
         */
        Real ref_chord;
        
        /*!
         *   returns the flight dynamic pressure
         */
        Real q0() const
        {
            return 0.5 * gas_property.rho * pow(velocity_magnitude(), 2);
        }
        
        Real p0()
        {
            return gas_property.pressure;
        }
        
        Real rho() const
        {
            return gas_property.rho;
        }
        
        /*!
         *    initializes the data structures
         */
        void init();
        
        Real rho_u1() const
        {
            return gas_property.rho * velocity_magnitude() * flow_unit_vector(0);
        }
        
        Real rho_u2() const
        {
            return gas_property.rho * velocity_magnitude() * flow_unit_vector(1);
        }
        
        Real rho_u3() const
        {
            return gas_property.rho * velocity_magnitude() * flow_unit_vector(2);
        }
        
        Real rho_e() const
        {
            return gas_property.rho * gas_property.cv * gas_property.T + q0();
        }

        Real rho_sens_rho() const {
            
            return 1.;
        }

        Real rho_u1_sens_rho() const {
            
            return this->rho_u1() / this->rho();
        }
        
        Real rho_u2_sens_rho() const {
            
            return this->rho_u2() / this->rho();;
        }
        
        Real rho_e_sens_rho() const
        {
            return this->rho_e() / this->rho();;
        }

        Real rho_u1_sens_mach() const {
            
            return this->rho_u1()/mach;
        }
        
        Real rho_u2_sens_mach() const {
            
            return this->rho_u2()/mach;
        }
        
        Real rho_e_sens_mach() const
        {
            return 2.*q0()/mach;
        }

        std::unique_ptr<MAST::FieldFunction<RealVectorX>> inf_sol;
    };
    
    
    
    inline
    void
    FlightCondition::init()
    {
        gas_property.init();
        flow_unit_vector /= flow_unit_vector.norm();
    }
    
}

#endif // __mast_fluid_condition_h__

