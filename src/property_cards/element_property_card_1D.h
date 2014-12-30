/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
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

#ifndef __mast__element_property_card_1D__
#define __mast__element_property_card_1D__


// MAST includes
#include "property_cards/element_property_card_base.h"



namespace MAST
{
    class ElementPropertyCard1D: public MAST::ElementPropertyCardBase {
        
    public:
        ElementPropertyCard1D():
        MAST::ElementPropertyCardBase(),
        _bending_model(MAST::DEFAULT_BENDING)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~ElementPropertyCard1D() { }
        
        
        /*!
         *   dimension of the element for which this property is defined
         */
        virtual unsigned int dim() const {
            return 1;
        }

        
        /*!
         *   returns the bending model to be used for the 2D element
         */
        void set_bending_model(MAST::BendingOperatorType b)  {
            _bending_model = b;
        }
        
        
        /*!
         *   returns the bending model to be used for the 2D element.
         */
        virtual MAST::BendingOperatorType bending_model(const libMesh::Elem& elem,
                                                        const libMesh::FEType& fe) const;
        
        
        /*!
         *    returns the extra quadrature order (on top of the system) that
         *    this element should use. This is elevated by two orders for a DKT
         *    element
         */
        virtual int extra_quadrature_order(const libMesh::Elem& elem,
                                           const libMesh::FEType& fe) const {
            if (this->bending_model(elem, fe) == MAST::BERNOULLI)
                return 2;
            else
                return 0;
        }
        
        
        /*!
         *   returns value of the property \par val. The string values for 
         *   \par val are IYY, IZZ, IYZ
         */
        virtual Real value(const std::string& val) const = 0;
        
        /*!
         *   vector in the x-y plane of the element. This should not be the same
         *   as the element x-axis.
         */
        libMesh::Point& y_vector() {
            return _local_y;
        }
        

        /*!
         *   constant reference to vector in the x-y plane of the element. 
         *   This should not be the same as the element x-axis.
         */
        const libMesh::Point& y_vector() const {
            return _local_y;
        }
        
        
    protected:
        
        /*!
         *   material property card. By default this chooses DKT for 3 noded
         *   triangles and Mindling for all other elements
         */
        MAST::BendingOperatorType _bending_model;
        
        /*!
         *   vector in the x-y plane.
         */
        libMesh::Point _local_y;
        
    };
    
    
}


#endif  // __mast__element_property_card_1D__
