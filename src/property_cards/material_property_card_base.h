/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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

#ifndef __mast__material_property_card_base__
#define __mast__material_property_card_base__

// MAST includes
#include "base/function_set_base.h"


namespace MAST
{
    // Forward decleration
    template <typename ValType> class FieldFunction;

    
    class  MaterialPropertyCardBase:
    public MAST::FunctionSetBase {
        
    public:
        
        MaterialPropertyCardBase():
        MAST::FunctionSetBase ()
        { }
        
        virtual ~MaterialPropertyCardBase() {}
        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        stiffness_matrix(const unsigned int dim,
                         const bool plane_stress = true) = 0;

        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        damping_matrix(const unsigned int dim) = 0;
        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        inertia_matrix(const unsigned int dim) = 0;

        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        thermal_expansion_matrix(const unsigned int dim) = 0;

        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        transverse_shear_stiffness_matrix() = 0;
        
        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        capacitance_matrix(const unsigned int dim) = 0;

        virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
        conductance_matrix(const unsigned int dim) = 0;

    protected:
        
    };
    
    
}


#endif // __mast__material_property_card_base__
