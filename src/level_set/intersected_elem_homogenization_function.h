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

#ifndef __mast__intersected_elem_homogenized_density_function_h__
#define __mast__intersected_elem_homogenized_density_function_h__

// MAST includes
#include "level_set/homogenized_density_function_base.h"


namespace MAST {

    // Forward declerations

    class LevelSetIntersection;

    class IntersectedElemHomogenizedDensityFunction:
    public MAST::HomogenizedDensityFunctionBase {
        
    public:
        
        IntersectedElemHomogenizedDensityFunction(const std::string& nm);
        
        virtual ~IntersectedElemHomogenizedDensityFunction();
        
        /*!
         *    computes and stores the volume fraction of each local element
         */
        virtual void initialize_element_volume_fractions();
        
        virtual void initialize_element_volume_fraction_sensitivity();
        
    protected:
        
        MAST::LevelSetIntersection           *_intersection;
        
    };
}

#endif // __mast__intersected_elem_homogenized_density_function_h__

