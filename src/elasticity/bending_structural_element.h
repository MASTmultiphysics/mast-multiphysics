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


#ifndef __mast__bending_structural_element__
#define __mast__bending_structural_element__

// MAST includes
#include "elasticity/structural_element_base.h"


namespace MAST {
    
    // Forward declerations
    class BoundaryConditionBase;
    class FEMOperatorMatrix;
    class GeomElem;
    
    
    class BendingStructuralElem:
    public MAST::StructuralElementBase {
        
    public:
        BendingStructuralElem(MAST::SystemInitialization& sys,
                              const MAST::GeomElem& elem,
                              const MAST::ElementPropertyCardBase& p);
        
        /*!
         *   destructor
         */
        virtual ~BendingStructuralElem();
        
        
        /*!
         *    row dimension of the direct strain matrix, also used for the
         *    bending operator row dimension
         */
        virtual unsigned int n_direct_strain_components() = 0;
        
        /*!
         *    row dimension of the von Karman strain matrix
         */
        virtual unsigned int n_von_karman_strain_components() = 0;
        
    protected:
        
    };
}


#endif // __mast__bending_structural_element__
