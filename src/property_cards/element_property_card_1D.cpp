
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

