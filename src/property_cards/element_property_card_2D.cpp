
// MAST includes
#include "property_cards/element_property_card_2D.h"



MAST::BendingOperatorType
MAST::ElementPropertyCard2D::bending_model(const libMesh::Elem& elem,
                                           const libMesh::FEType& fe) const {
    // for a TRI3 element, default bending is DKT. For all other elements
    // the default is Mindlin. Otherwise it returns the model set for
    // this card.
    switch (elem.type()) {
        case libMesh::TRI3:
            if ((fe.family == libMesh::LAGRANGE) &&
                (fe.order  == libMesh::FIRST) &&
                (_bending_model == MAST::DEFAULT_BENDING))
                return MAST::DKT;
            else
                return MAST::MINDLIN;
            break;
            
        default:
            if (_bending_model == MAST::DEFAULT_BENDING)
                return MAST::MINDLIN;
            else
                return _bending_model;
            break;
    }
}

