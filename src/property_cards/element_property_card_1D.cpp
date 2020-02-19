/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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
#include "mesh/geom_elem.h"

MAST::BendingOperatorType
MAST::ElementPropertyCard1D::bending_model(const MAST::GeomElem& elem) const {
    
    // Ensure a valid bending_model is chosen, raise error if not
    if ((_bending_model == MAST::DKT) or (_bending_model == MAST::MINDLIN)){
        libmesh_error_msg("Invalid bending model for 1D element. Should be either MAST::BERNOULLI, MAST::TIMOSHENKO, MAST::NO_BENDING, OR MAST::DEFAULT_BENDING; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
    }
    
    // for an EDGE2 element, default bending is Bernoulli. For all other elements
    // the default is Timoshenko. Otherwise it returns the model set for
    // this card.
    switch (elem.get_reference_elem().type()) {
        case libMesh::EDGE2:
            // assuming that all variables have the same interpolation
            if ((elem.get_fe_type(0).family == libMesh::LAGRANGE) &&
                (elem.get_fe_type(0).order  == libMesh::FIRST) &&
                (_bending_model == MAST::DEFAULT_BENDING))
                return MAST::BERNOULLI;
            else
                return _bending_model;
            break;
            
        default:
            if (_bending_model == MAST::DEFAULT_BENDING)
                return MAST::TIMOSHENKO;
            else
                return _bending_model;
            break;
    }
}

