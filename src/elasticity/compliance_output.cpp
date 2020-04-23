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
#include "elasticity/compliance_output.h"
#include "elasticity/structural_element_base.h"
#include "base/assembly_base.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "base/physics_discipline_base.h"
#include "base/boundary_condition_base.h"
#include "property_cards/element_property_card_1D.h"
#include "level_set/level_set_intersection.h"
#include "level_set/level_set_intersected_elem.h"
#include "mesh/geom_elem.h"

// libMesh includes
#include "libmesh/parallel.h"

MAST::ComplianceOutput::ComplianceOutput():
MAST::OutputAssemblyElemOperations(),
_compliance       (0.),
_dcompliance_dp   (0.)  {
    
}




MAST::ComplianceOutput::~ComplianceOutput() {

}




void
MAST::ComplianceOutput::zero_for_analysis() {
    
    _compliance         = 0.;
    _dcompliance_dp     = 0.;
}



void
MAST::ComplianceOutput::zero_for_sensitivity() {
    
    _compliance         = 0.;
    _dcompliance_dp     = 0.;
}


void
MAST::ComplianceOutput::evaluate() {
    
    // make sure that this has not been initialized ana calculated for all elems
    libmesh_assert(_physics_elem);
    
    if (this->if_evaluate_for_element(_physics_elem->elem())) {
        
        MAST::StructuralElementBase& e =
        dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        
        RealVectorX
        vec   = RealVectorX::Zero(e.sol().size());
        
        RealMatrixX
        dummy = RealMatrixX::Zero(vec.size(), vec.size());
        
        e.side_external_residual(false,
                                 vec,
                                 dummy,
                                 dummy,
                                 _discipline->side_loads());
        e.volume_external_residual(false,
                                   vec,
                                   dummy,
                                   dummy,
                                   _discipline->volume_loads());

        // compute the contribution of this element to compliance
        _compliance -= vec.dot(e.sol());
    }
}



void
MAST::ComplianceOutput::evaluate_sensitivity(const MAST::FunctionBase &f) {
    
    // make sure that this has not been initialized ana calculated for all elems
    libmesh_assert(_physics_elem);
    
    if (this->if_evaluate_for_element(_physics_elem->elem())) {
        
        MAST::StructuralElementBase& e =
        dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        
        RealVectorX
        vec   = RealVectorX::Zero(e.sol().size());
        
        RealMatrixX
        dummy = RealMatrixX::Zero(vec.size(), vec.size());
        
        e.side_external_residual_sensitivity(f, false,
                                             vec,
                                             dummy,
                                             dummy,
                                             _discipline->side_loads());
        e.volume_external_residual_sensitivity(f, false,
                                               vec,
                                               dummy,
                                               dummy,
                                               _discipline->volume_loads());
        
        // ask for the values
        _dcompliance_dp -= vec.dot(e.sol());
        
        vec.setZero();
        e.side_external_residual(false,
                                 vec,
                                 dummy,
                                 dummy,
                                 _discipline->side_loads());
        e.volume_external_residual(false,
                                   vec,
                                   dummy,
                                   dummy,
                                   _discipline->volume_loads());
        
        // ask for the values
        _dcompliance_dp -= vec.dot(e.sol(true));
    }
}



void
MAST::ComplianceOutput::
evaluate_topology_sensitivity(const MAST::FunctionBase &f) {
    
    // the primal data should have been calculated
    libmesh_assert(_physics_elem);
    libmesh_assert(f.is_topology_parameter());
    
    if (this->if_evaluate_for_element(_physics_elem->elem())) {
        
        std::pair<const MAST::FieldFunction<RealVectorX>*, unsigned int>
        val = this->get_elem_boundary_velocity_data();
        
        if (val.first) {

            MAST::StructuralElementBase& e =
            dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
            
            RealVectorX
            vec   = RealVectorX::Zero(e.sol().size());
            
            RealMatrixX
            dummy = RealMatrixX::Zero(vec.size(), vec.size());
            
            e.volume_external_residual_boundary_velocity(f,
                                                         val.second,
                                                         *val.first,
                                                         _discipline->volume_loads(),
                                                         false,
                                                         vec,
                                                         dummy);
            
            // compute the contribution of this element to compliance
            _dcompliance_dp -= vec.dot(e.sol());
        }
    }
}



void
MAST::ComplianceOutput::
evaluate_topology_sensitivity(const MAST::FunctionBase &f,
                              const MAST::FieldFunction<RealVectorX> &vel) {
    
    // the primal data should have been calculated
    libmesh_assert(_physics_elem);
    libmesh_assert(f.is_topology_parameter());

    const MAST::LevelSetIntersectedElem
    &elem = dynamic_cast<const MAST::LevelSetIntersectedElem&>(_physics_elem->elem());
    
    // sensitivity only exists at the boundary. So, we proceed with calculation
    // only if this element has an intersection in the interior, or with a side.
    
    if (this->if_evaluate_for_element(elem)    &&
        elem.if_elem_has_level_set_boundary()       &&
        elem.if_subelem_has_side_on_level_set_boundary()) {
        
        MAST::StructuralElementBase& e =
        dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        
        RealVectorX
        vec   = RealVectorX::Zero(e.sol().size());
        
        RealMatrixX
        dummy = RealMatrixX::Zero(vec.size(), vec.size());
        
        e.volume_external_residual_boundary_velocity(f,
                                                     elem.get_subelem_side_on_level_set_boundary(),
                                                     vel,
                                                     _discipline->volume_loads(),
                                                     false,
                                                     vec,
                                                     dummy);
        
        // compute the contribution of this element to compliance
        _dcompliance_dp -= vec.dot(e.sol());
    }
}



Real
MAST::ComplianceOutput::output_total() {
    
    Real val = _compliance;
    
    if (!_skip_comm_sum)
        _system->system().comm().sum(val);
    
    return val;
}



Real
MAST::ComplianceOutput::output_sensitivity_total(const MAST::FunctionBase& p) {
    
    Real val = _dcompliance_dp;
    
    if (!_skip_comm_sum)
        _system->system().comm().sum(val);

    return val;
}



void
MAST::ComplianceOutput::output_derivative_for_elem(RealVectorX& dq_dX) {
    
    // make sure that this has not been initialized ana calculated for all elems
    libmesh_assert(_physics_elem);
    
     // since compliance = -X^T F,  derivative wrt X = -F.
    
    if (this->if_evaluate_for_element(_physics_elem->elem())) {
        
        MAST::StructuralElementBase& e =
        dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        
        dq_dX.setZero();
        
        RealMatrixX
        dummy = RealMatrixX::Zero(dq_dX.size(), dq_dX.size());
        
        e.side_external_residual(false,
                                 dq_dX,
                                 dummy,
                                 dummy,
                                 _discipline->side_loads());
        e.volume_external_residual(false,
                                   dq_dX,
                                   dummy,
                                   dummy,
                                   _discipline->volume_loads());
        
        dq_dX *= -1.;
    }
}



void
MAST::ComplianceOutput::
set_elem_data(unsigned int dim,
              const libMesh::Elem& ref_elem,
              MAST::GeomElem& elem) const {
    
    libmesh_assert(!_physics_elem);
    
    if (dim == 1) {
        
        const MAST::ElementPropertyCard1D& p =
        dynamic_cast<const MAST::ElementPropertyCard1D&>(_discipline->get_property_card(ref_elem));
        
        elem.set_local_y_vector(p.y_vector());
    }
}


void
MAST::ComplianceOutput::init(const MAST::GeomElem& elem) {
    
    libmesh_assert(!_physics_elem);
    libmesh_assert(_assembly);
    libmesh_assert(_system);
    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    _physics_elem =
    MAST::build_structural_element(*_system, elem, p).release();
}


