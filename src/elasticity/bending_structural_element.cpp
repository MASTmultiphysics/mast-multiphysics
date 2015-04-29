
// MAST includes
#include "elasticity/bending_structural_element.h"
#include "elasticity/bending_operator.h"
#include "property_cards/element_property_card_base.h"


MAST::BendingStructuralElem::
BendingStructuralElem(MAST::SystemInitialization& sys,
                      const libMesh::Elem& elem,
                      const MAST::ElementPropertyCardBase& p,
                      const bool output_eval_mode):
MAST::StructuralElementBase(sys, elem, p, output_eval_mode) {
    
    // initialize the bending operator
    MAST::BendingOperatorType bending_model =
    _property.bending_model(elem, _fe->get_fe_type());
    
    _bending_operator.reset(MAST::build_bending_operator(bending_model, *this).release());
}
