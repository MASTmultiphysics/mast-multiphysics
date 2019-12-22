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

#ifndef __mast__heaviside_elem_homogenized_density_function_h__
#define __mast__heaviside_elem_homogenized_density_function_h__

// MAST includes
#include "level_set/homogenized_density_function_base.h"


namespace MAST {
    
    class HeavisideElemHomogenizedDensityFunction:
    public MAST::HomogenizedDensityFunctionBase {
        
    public:
        
        HeavisideElemHomogenizedDensityFunction(const std::string& nm);
        
        virtual ~HeavisideElemHomogenizedDensityFunction();
        
        /*!
         *   width over which the approximate Heaviside function is smoothed. Note that this is the
         *   width fo the level-set function value, and not necessarily the width of the geometric element.
         *   For level-set functions that are initialized to have a unit gradient norm, these two would
         *   coincide. The default value is 0.01.
         */
        void set_smoothing_width(Real w) { _width = w;}
        
        virtual void initialize_element_volume_fractions();
        
        virtual void initialize_element_volume_fraction_sensitivity();
        
    protected:
        
        
        Real  _width;
    };
}

#endif // __mast__heaviside_elem_homogenized_density_function_h__

