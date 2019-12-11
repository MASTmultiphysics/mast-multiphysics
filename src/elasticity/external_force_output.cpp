/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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
#include "elasticity/external_force_output.h"
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
#include "libmesh/dof_map.h"


MAST::ExternalForceOutput::ExternalForceOutput(const std::vector<Real> w):
MAST::OutputAssemblyElemOperations(),
_external_force       (0.),
_dexternal_force_dp   (0.),
_w                    (w){
}




MAST::ExternalForceOutput::~ExternalForceOutput() {

}




void
MAST::ExternalForceOutput::zero_for_analysis() {
    
    _external_force         = 0.;
    _dexternal_force_dp     = 0.;
}



void
MAST::ExternalForceOutput::zero_for_sensitivity() {
    
    _external_force         = 0.;
    _dexternal_force_dp     = 0.;
}


void MAST::ExternalForceOutput::evaluate() 
{
    // make sure that this has not been initialized ana calculated for all elems
    libmesh_assert(_physics_elem);
    
    if (this->if_evaluate_for_element(_physics_elem->elem())) 
    {
        MAST::StructuralElementBase& e =
            dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        
        RealVectorX vec   = RealVectorX::Zero(e.sol().size());
        
        RealMatrixX dummy = RealMatrixX::Zero(vec.size(), vec.size());
        
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
        
        const libMesh::DofMap& dof_map = _system->system().get_dof_map();
        std::vector<libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(&_physics_elem->elem().get_reference_elem(), 
                            dof_indices);
        
        RealVectorX w = RealVectorX::Zero(dof_indices.size());
        // compute the contribution of this element to internal_force combo
        for (uint i=0; i<dof_indices.size(); i++)
        {
            w(i) = _w[dof_indices[i]];
        }

        // compute the contribution of this element to external_force
        _external_force -= vec.dot(w);
    }
}


void MAST::ExternalForceOutput::evaluate_for_node(const RealVectorX& Xnode,
                                               const RealVectorX& Fpnode)
{    
    // make sure that this has not been initialized and calculated for all nodes
    libmesh_assert(_node);
    
    if (this->if_evaluate_for_node(*_node)) 
    {
        const libMesh::DofMap& dof_map = _system->system().get_dof_map();
        std::vector<libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(_node, dof_indices);
         
        RealVectorX w = RealVectorX::Zero(dof_indices.size());
        for (uint i=0; i<Xnode.size(); i++)
        {
            w(i) = _w[dof_indices[i]];
        }
        
        _external_force += w.dot(Fpnode);
    }
}


void MAST::ExternalForceOutput::evaluate_sensitivity(const MAST::FunctionBase &f) 
{
    // make sure that this has not been initialized ana calculated for all elems
    libmesh_assert(_physics_elem);
    
    if (this->if_evaluate_for_element(_physics_elem->elem())) 
    {
        MAST::StructuralElementBase& e =
            dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        
        RealVectorX vec   = RealVectorX::Zero(e.sol().size());
        
        RealMatrixX dummy = RealMatrixX::Zero(vec.size(), vec.size());
        
        const libMesh::DofMap& dof_map = _system->system().get_dof_map();
        std::vector<libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(&_physics_elem->elem().get_reference_elem(), 
                            dof_indices);
        
        RealVectorX w = RealVectorX::Zero(dof_indices.size());
        // compute the contribution of this element to internal_force combo
        for (uint i=0; i<dof_indices.size(); i++)
        {
            w(i) = _w[dof_indices[i]];
        }
        
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
        _dexternal_force_dp -= vec.dot(w);
        
        RealVectorX dq_dX = RealVectorX::Zero(vec.size());
        output_derivative_for_elem(dq_dX);
        
        // ask for the values
        _dexternal_force_dp -= dq_dX.dot(e.sol(true));
    }
}


void MAST::ExternalForceOutput::evaluate_sensitivity_for_node(
    const MAST::FunctionBase& f, const RealVectorX& Xnode, 
    const RealVectorX& dXnode_dparam, const RealVectorX& Fpnode, 
    const RealVectorX& dpF_fpparam_node)
{
    // make sure that this has not been initialized and calculated for all nodes
    libmesh_assert(_node);
    
    if (this->if_evaluate_for_node(*_node)) 
    {
        RealVectorX dpr_dpXnode = RealVectorX::Zero(Xnode.size());
        output_derivative_for_node(Xnode, Fpnode, dpr_dpXnode);
        
        _dexternal_force_dp += dpF_fpparam_node.dot(Xnode) + dpr_dpXnode.dot(dXnode_dparam);
    }
}


void MAST::ExternalForceOutput::evaluate_topology_sensitivity(
    const MAST::FunctionBase &f, const MAST::FieldFunction<RealVectorX> &vel) 
{
    libmesh_error_msg("evaluate_topology_sensitivity is not implemeneted in ExternalForceOutput");
}



Real
MAST::ExternalForceOutput::output_total() {
    
    Real val = _external_force;
    _system->system().comm().sum(val);
    return val;
}



Real MAST::ExternalForceOutput::output_sensitivity_total(const MAST::FunctionBase& p) 
{
    Real val = _dexternal_force_dp;
    _system->system().comm().sum(val);
    return val;
}



void MAST::ExternalForceOutput::output_derivative_for_elem(RealVectorX& dq_dX) 
{
    // make sure that this has not been initialized and calculated for all elems
    libmesh_assert(_physics_elem);
    
    if (this->if_evaluate_for_element(_physics_elem->elem())) 
    {
        MAST::StructuralElementBase& e =
            dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        
        dq_dX.setZero();
        
        RealMatrixX dummy = RealMatrixX::Zero(dq_dX.size(), dq_dX.size());
        
        RealMatrixX dpfext_dpX = RealMatrixX::Zero(dq_dX.size(), dq_dX.size());
        
        e.side_external_residual(true,
                                 dq_dX,
                                 dummy,
                                 dpfext_dpX,
                                 _discipline->side_loads());
        e.volume_external_residual(true,
                                   dq_dX,
                                   dummy,
                                   dpfext_dpX,
                                   _discipline->volume_loads());
        
        const libMesh::DofMap& dof_map = _system->system().get_dof_map();
        std::vector<libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(&_physics_elem->elem().get_reference_elem(), 
                            dof_indices);
        
        RealVectorX w = RealVectorX::Zero(dof_indices.size());
        // compute the contribution of this element to internal_force combo
        for (uint i=0; i<dof_indices.size(); i++)
        {
            w(i) = _w[dof_indices[i]];
        }
        
        dq_dX += dpfext_dpX.transpose() * w;
    }
}


void MAST::ExternalForceOutput::output_derivative_for_node(
    const RealVectorX& Xnode, const RealVectorX& Fpnode, RealVectorX& dq_dX)
{
    // make sure that this has not been initialized and calculated for all nodes
    libmesh_assert(_node);
    
    if (this->if_evaluate_for_node(*_node)) 
    {
        dq_dX.setZero();
        // TODO: Update this if follower point loads are implemented where point loads will depend on displacement.
    }
}


void
MAST::ExternalForceOutput::
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
MAST::ExternalForceOutput::init(const MAST::GeomElem& elem) {
    
    libmesh_assert(!_physics_elem);
    libmesh_assert(_assembly);
    libmesh_assert(_system);
    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    _physics_elem =
    MAST::build_structural_element(*_system, *_assembly, elem, p).release();
}


