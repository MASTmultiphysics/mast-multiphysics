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

#ifndef __mast__isotropic_material_property_card__
#define __mast__isotropic_material_property_card__

// MAST includes
#include "property_cards/material_property_card_base.h"


namespace MAST {
    

    class IsotropicMaterialPropertyCard:
    public MAST::MaterialPropertyCardBase {
        
    public:
        
        IsotropicMaterialPropertyCard();
        
        virtual ~IsotropicMaterialPropertyCard();
        
        
        virtual const MAST::FieldFunction<RealMatrixX>&
        stiffness_matrix(const unsigned int dim,
                         const bool plane_stress = true);
        
        virtual const MAST::FieldFunction<RealMatrixX>&
        damping_matrix(const unsigned int dim);
        
        virtual const MAST::FieldFunction<RealMatrixX>&
        inertia_matrix(const unsigned int dim);
        
        virtual const MAST::FieldFunction<RealMatrixX>&
        thermal_expansion_matrix(const unsigned int dim);
        
        virtual const MAST::FieldFunction<RealMatrixX>&
        transverse_shear_stiffness_matrix();
        
        virtual const MAST::FieldFunction<RealMatrixX>&
        capacitance_matrix(const unsigned int dim);
        
        virtual const MAST::FieldFunction<RealMatrixX>&
        conductance_matrix(const unsigned int dim);
        
    protected:
        
        MAST::FieldFunction<RealMatrixX>
        *_stiff_mat_1d,
        *_stiff_mat_2d,
        *_stiff_mat_3d,
        *_damp_mat,
        *_inertia_mat_3d,
        *_thermal_exp_mat_1d,
        *_thermal_exp_mat_2d,
        *_thermal_exp_mat_3d,
        *_transverse_shear_mat,
        *_thermal_capacitance_mat_1d,
        *_thermal_capacitance_mat_2d,
        *_thermal_capacitance_mat_3d,
        *_thermal_conductance_mat_1d,
        *_thermal_conductance_mat_2d,
        *_thermal_conductance_mat_3d;
    };
    
}





#endif // __mast__isotropic_material_property_card__
