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

// MAST includes
#include "property_cards/element_property_card_1D.h"


MAST::BendingOperatorType
MAST::ElementPropertyCard1D::bending_model(const libMesh::Elem& elem,
                                           const libMesh::FEType& fe) const {
    // for an EDGE2 element, default bending is Bernoulli. For all other elements
    // the default is Timoshenko. Otherwise it returns the model set for
    // this card.
    switch (elem.type()) {
        case libMesh::EDGE2:
            if ((fe.family == libMesh::LAGRANGE) &&
                (fe.order  == libMesh::FIRST) &&
                (_bending_model == MAST::DEFAULT_BENDING))
                return MAST::BERNOULLI;
            else
                return MAST::TIMOSHENKO;
            break;
            
        default:
            if (_bending_model == MAST::DEFAULT_BENDING)
                return MAST::TIMOSHENKO;
            else
                return _bending_model;
            break;
    }
}

