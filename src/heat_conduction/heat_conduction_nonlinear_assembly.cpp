/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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
#include "heat_conduction/heat_conduction_nonlinear_assembly.h"
#include "base/assembly_base.h"
#include "heat_conduction/heat_conduction_elem_base.h"
#include "property_cards/element_property_card_1D.h"
#include "base/physics_discipline_base.h"
#include "mesh/local_elem_fe.h"


MAST::HeatConductionNonlinearAssemblyElemOperations::
HeatConductionNonlinearAssemblyElemOperations():
MAST::NonlinearImplicitAssemblyElemOperations() {
    
}



MAST::HeatConductionNonlinearAssemblyElemOperations::
~HeatConductionNonlinearAssemblyElemOperations() {
    
}



std::unique_ptr<MAST::ElementBase>
MAST::HeatConductionNonlinearAssemblyElemOperations::
build_elem(const libMesh::Elem& elem) {

    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>
    (_assembly->discipline().get_property_card(elem));
    
    MAST::ElementBase* rval =
    new MAST::HeatConductionElementBase(_assembly->system_init(), *_assembly, elem, p);
    
    return std::unique_ptr<MAST::ElementBase>(rval);
}




void
MAST::HeatConductionNonlinearAssemblyElemOperations::
elem_calculations(MAST::ElementBase& elem,
                  bool if_jac,
                  RealVectorX& vec,
                  RealMatrixX& mat) {
    
    MAST::HeatConductionElementBase& e =
    dynamic_cast<MAST::HeatConductionElementBase&>(elem);
    
    vec.setZero();
    mat.setZero();
    
    e.internal_residual(if_jac, vec, mat);
    e.side_external_residual(if_jac, vec, mat, _assembly->discipline().side_loads());
    e.volume_external_residual(if_jac, vec, mat, _assembly->discipline().volume_loads());
}




void
MAST::HeatConductionNonlinearAssemblyElemOperations::
elem_sensitivity_calculations(MAST::ElementBase& elem,
                              RealVectorX& vec) {
    
    MAST::HeatConductionElementBase& e =
    dynamic_cast<MAST::HeatConductionElementBase&>(elem);
    
    vec.setZero();
    RealMatrixX
    dummy = RealMatrixX::Zero(vec.size(), vec.size());
    
    e.internal_residual_sensitivity(false, vec, dummy);
    e.side_external_residual_sensitivity(false, vec, dummy, _assembly->discipline().side_loads());
    e.volume_external_residual_sensitivity(false, vec, dummy, _assembly->discipline().volume_loads());
}



void
MAST::HeatConductionNonlinearAssemblyElemOperations::
elem_second_derivative_dot_solution_assembly(MAST::ElementBase& elem,
                                             RealMatrixX& m) {
    
    libmesh_error(); // to be implemented
}


void
MAST::HeatConductionNonlinearAssemblyElemOperations::
set_local_fe_data(const libMesh::Elem& e,
                  MAST::LocalElemFE& fe) const {
    
    if (e.dim() == 1) {
        
        const MAST::ElementPropertyCard1D&
        p_card = dynamic_cast<const MAST::ElementPropertyCard1D&>
        (_assembly->discipline().get_property_card(e));
        
        fe.set_1d_y_vector(p_card.y_vector());
    }
}
