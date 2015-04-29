

// MAST includes
#include "elasticity/structural_transient_assembly.h"
#include "elasticity/structural_element_base.h"
#include "property_cards/element_property_card_base.h"
#include "base/physics_discipline_base.h"


MAST::StructuralTransientAssembly::
StructuralTransientAssembly():
MAST::TransientAssembly() {
    
}




MAST::StructuralTransientAssembly::
~StructuralTransientAssembly() {
    
}



void
MAST::StructuralTransientAssembly::
_elem_calculations(MAST::ElementBase& elem,
                   bool if_jac,
                   RealVectorX& f_m,
                   RealVectorX& f_x,
                   RealMatrixX& f_m_jac_xdot,
                   RealMatrixX& f_m_jac,
                   RealMatrixX& f_x_jac) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    f_m.setZero();
    f_x.setZero();
    f_m_jac_xdot.setZero();
    f_m_jac.setZero();
    f_x_jac.setZero();
    
    // assembly of the flux terms
    e.internal_residual(if_jac, f_x, f_x_jac, false);
    e.side_external_residual<Real>(if_jac, f_x, f_x_jac, _discipline->side_loads());
    e.volume_external_residual<Real>(if_jac, f_x, f_x_jac, _discipline->volume_loads());
    
    //assembly of the capacitance term
    e.damping_residual(if_jac, f_m, f_m_jac_xdot, f_m_jac);
}




void
MAST::StructuralTransientAssembly::
_elem_calculations(MAST::ElementBase& elem,
                   bool if_jac,
                   RealVectorX& f_m,
                   RealVectorX& f_x,
                   RealMatrixX& f_m_jac_xddot,
                   RealMatrixX& f_m_jac_xdot,
                   RealMatrixX& f_m_jac,
                   RealMatrixX& f_x_jac_xdot,
                   RealMatrixX& f_x_jac) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    f_m.setZero();
    f_x.setZero();
    f_m_jac_xddot.setZero();
    f_m_jac_xdot.setZero();
    f_m_jac.setZero();
    f_x_jac_xdot.setZero();
    f_x_jac.setZero();
    
    // assembly of the flux terms
    e.internal_residual(if_jac, f_x, f_x_jac, false);
    e.damping_residual(if_jac, f_x, f_x_jac_xdot, f_x_jac);
    e.side_external_residual<Real>(if_jac, f_x, f_x_jac, _discipline->side_loads());
    e.volume_external_residual<Real>(if_jac, f_x, f_x_jac, _discipline->volume_loads());
    
    //assembly of the capacitance term
    e.inertial_residual(if_jac, f_m, f_m_jac_xddot, f_m_jac_xdot, f_m_jac);
}



void
MAST::StructuralTransientAssembly::
_elem_sensitivity_calculations(MAST::ElementBase& elem,
                               RealVectorX& vec) {
    
    
    
}



std::auto_ptr<MAST::ElementBase>
MAST::StructuralTransientAssembly::_build_elem(const libMesh::Elem& elem) {
    
    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    MAST::ElementBase* rval =
    MAST::build_structural_element(*_system, elem, p, false).release();
    
    return std::auto_ptr<MAST::ElementBase>(rval);
}
