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

#ifndef __mast__element_property_card_2D__
#define __mast__element_property_card_2D__

// MAST includes
#include "property_cards/element_property_card_base.h"


namespace MAST
{
    
    class ElementPropertyCard2D: public MAST::ElementPropertyCardBase {
        
    public:
        ElementPropertyCard2D():
        MAST::ElementPropertyCardBase(),
        _bending_model(MAST::DEFAULT_BENDING),
        _if_plane_stress(true)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~ElementPropertyCard2D() { }
        
        /*!
         *   returns the bending model to be used for the 2D element
         */
        void set_bending_model(MAST::BendingOperatorType b)  {
            _bending_model = b;
        }
        
        
        /*!
         *   returns the bending model to be used for the 2D element.
         */
        MAST::BendingOperatorType bending_model(const libMesh::Elem& elem,
                                         const libMesh::FEType& fe) const;
        
        
        /*!
         *    returns the extra quadrature order (on top of the system) that
         *    this element should use. This is elevated by two orders for a DKT
         *    element
         */
        virtual int extra_quadrature_order(const libMesh::Elem& elem,
                                           const libMesh::FEType& fe) const {
            if (this->bending_model(elem, fe) == MAST::DKT)
                return 2;
            else
                return 0;
        }
        
        /*!
         *   sets the flag for plane stress.
         */
        void set_plane_stress(bool val) {
            _if_plane_stress = val;
        }
        
        /*!
         *   returns the flag for plane stress.
         */
        bool plane_stress() const {
            return _if_plane_stress;
        }
        
        
    protected:
        
        /*!
         *   material property card. By default this chooses DKT for 3 noded
         *   triangles and Mindling for all other elements
         */
        MAST::BendingOperatorType _bending_model;
        
        /*!
         *   if the analysis is plne stress, otherwise it is plane strain.
         *   Note that this is true by default
         */
        bool _if_plane_stress;
        
    };
    
}




#endif // __mast__element_property_card_2D__
