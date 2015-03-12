

// MAST includes
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_element_base.h"
#include "property_cards/element_property_card_base.h"
#include "base/physics_discipline_base.h"



MAST::StructuralNonlinearAssembly::
StructuralNonlinearAssembly():
MAST::NonlinearImplicitAssembly() {
    
}



MAST::StructuralNonlinearAssembly::
~StructuralNonlinearAssembly() {
    
}



std::auto_ptr<MAST::ElementBase>
MAST::StructuralNonlinearAssembly::_build_elem(const libMesh::Elem& elem) {
    
    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    MAST::ElementBase* rval =
    MAST::build_structural_element(*_system, elem, p).release();
    
    return std::auto_ptr<MAST::ElementBase>(rval);
}




void
MAST::StructuralNonlinearAssembly::
_elem_calculations(MAST::ElementBase& elem,
                   bool if_jac,
                   RealVectorX& vec,
                   RealMatrixX& mat) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    vec.setZero();
    mat.setZero();
    
    e.internal_residual(if_jac, vec, mat, false);
    e.side_external_residual<Real>(if_jac, vec, mat, _discipline->side_loads());
    e.volume_external_residual<Real>(if_jac, vec, mat, _discipline->volume_loads());
}




void
MAST::StructuralNonlinearAssembly::
_elem_sensitivity_calculations(MAST::ElementBase& elem,
                               RealVectorX& vec) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    vec.setZero();
    RealMatrixX mat; // dummy matrix

    libmesh_error();
//    e.internal_residual_sensitivity(false, vec, mat, false);
//    e.side_external_residual_sensitivity<Real>(false, vec, mat, _discipline->side_loads());
//    e.volume_external_residual_sensitivity<Real>(false, vec, mat, _discipline->volume_loads());
}





