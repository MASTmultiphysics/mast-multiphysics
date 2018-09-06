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

#ifndef __mast__frequency_function_h__
#define __mast__frequency_function_h__

// MAST includes
#include "base/field_function_base.h"


namespace MAST {
    
    // Forward declerations
    class Parameter;
    
    class FrequencyFunction:
    public MAST::FieldFunction<Real> {
    
    public:
        
        FrequencyFunction(const std::string& nm,
                          MAST::FieldFunction<Real>& omega,
                          MAST::FieldFunction<Real>& velocity,
                          MAST::FieldFunction<Real>& b_ref);
        
        
        virtual ~FrequencyFunction();
        
        
        /*!
         *   @returns the flag for whether or not the frequency has been
         *   nondimensionalized: omega b/V
         */
        bool if_nondimensional() const {
            
            return _if_red_freq;
        }

        
        /*!
         *   sets the flag for whether or not the frequency has been
         *   nondimensionalized: omega b/V
         */
        inline void if_nondimensional(bool v)  {
            
            _if_red_freq = v;
        }

        
        /*!
         *    calculates the value of the function and returns it in \p v.
         */
        virtual void operator() (Real& v) const;
        
        
        /*!
         *    calculates the value of the function derivative and
         *    returns it in \p v.
         */
        virtual void derivative (const MAST::FunctionBase& f,
                                 Real& v) const;

        
        /*!
         *    @returns value of the non-dimensionalizing factor (b/V) to
         *    convert frequency (omega) to reduced frequency (omega * b/V)
         */
        void nondimensionalizing_factor(Real& v);
         
        
        
    protected:

        
        bool _if_red_freq;
        
        MAST::FieldFunction<Real>
        &_omega,
        &_velocity,
        &_b_ref;
        
    };
}

#endif // __mast__frequency_function_h__
