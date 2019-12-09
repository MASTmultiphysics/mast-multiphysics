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

#include <fstream>

// MAST includes
#include "elasticity/displacement_output.h"
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
#include "libmesh/dof_map.h"


MAST::DisplacementOutput::DisplacementOutput(const std::vector<Real> w):
    MAST::OutputAssemblyElemOperations(),
    _displacement       (0.), _ddisplacement_dp   (0.), _w(w)
{
}


MAST::DisplacementOutput::~DisplacementOutput() 
{
}


void MAST::DisplacementOutput::zero_for_analysis() 
{
    
    _displacement         = 0.;
    _ddisplacement_dp     = 0.;
}


void MAST::DisplacementOutput::zero_for_sensitivity() 
{
    
    _displacement         = 0.;
    _ddisplacement_dp     = 0.;
}


/**
 * This is called by the calculate_output method in the assembly object.
 */
void MAST::DisplacementOutput::evaluate() 
{
    // make sure that this has not been initialized and calculated for all elems
    libmesh_assert(_physics_elem);
    
    // This output is node based only. It is not an element quantity.
}


void MAST::DisplacementOutput::evaluate_for_node(const RealVectorX& Xnode,
                                                 const RealVectorX& Fpnode)
{    
    // make sure that this has not been initialized and calculated for all nodes
    libmesh_assert(_node);
    
    if (this->if_evaluate_for_node(*_node)) 
    {
        const libMesh::DofMap& dof_map = _system->system().get_dof_map();
        std::vector<libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(_node, dof_indices);
                
        for (uint i=0; i<Xnode.size(); i++)
        {
            _displacement += _w[dof_indices[i]] * Xnode(i);
        }
    }
}


/**
 * This is used by the calculate_output_direct_sensitivity method in the 
 * assembly object. It calculates dpResponseElement_dpparameter.
 */
void MAST::DisplacementOutput::evaluate_sensitivity(const MAST::FunctionBase &f)
{
    
    // make sure that this has not been initialized ana calculated for all elems
    libmesh_assert(_physics_elem);
    
    // Nothing to do here, _ddisplacement_dp is assumed to always be zero
    // (since displacements shouldn't normally explicitly depend on parameters)
    // and _ddisplacement_dp is already set to zero initially.
}


void MAST::DisplacementOutput::evaluate_sensitivity_for_node(
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
        
        _ddisplacement_dp += dpr_dpXnode.dot(dXnode_dparam);
    }
}


Real MAST::DisplacementOutput::output_total() 
{
    Real val = _displacement;
    _system->system().comm().sum(val);
    return val;
}


/**
 * This is calculates dpRepsonse_dpParameter and is used by the
 * calculate_output_adjoint sensitivity method in the assembly object.
 */
Real MAST::DisplacementOutput::output_sensitivity_total(
    const MAST::FunctionBase& p) 
{
    Real val = _ddisplacement_dp;
    _system->system().comm().sum(val);
    return val;
}


/**
 * This is calculating dpResponse_dpStateVector for a single element.
 * This gets called by calculate_output_derivative in the assembly object
 * to assemble the total sensitivity.
 */
void MAST::DisplacementOutput::output_derivative_for_elem(RealVectorX& dq_dX) 
{
    // make sure that this has not been initialized and calculated for all elems
    libmesh_assert(_physics_elem);
        
    // Nothing to do here, displacement is a nodal quantitity and is    
    // not an element quantity
}


void MAST::DisplacementOutput::output_derivative_for_node(
    const RealVectorX& Xnode, const RealVectorX& Fpnode, RealVectorX& dq_dX)
{
    // make sure that this has not been initialized and calculated for all nodes
    libmesh_assert(_node);
    
    if (this->if_evaluate_for_node(*_node)) 
    {
        dq_dX.setZero();
        
        const libMesh::DofMap& dof_map = _system->system().get_dof_map();
        std::vector<libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(_node, dof_indices);
        
        for (uint i=0; i<dof_indices.size(); i++)
        {
            dq_dX(i) = _w[dof_indices[i]];
        }
    }
}


void MAST::DisplacementOutput:: set_elem_data(unsigned int dim,
    const libMesh::Elem& ref_elem, MAST::GeomElem& elem) const 
{
    libmesh_assert(!_physics_elem);
    
    if (dim == 1) 
    {
        const MAST::ElementPropertyCard1D& p =
        dynamic_cast<const MAST::ElementPropertyCard1D&>(_discipline->get_property_card(ref_elem));
        
        elem.set_local_y_vector(p.y_vector());
    }
}


void MAST::DisplacementOutput::init(const MAST::GeomElem& elem) 
{
    libmesh_assert(!_physics_elem);
    libmesh_assert(_assembly);
    libmesh_assert(_system);
    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    _physics_elem =
    MAST::build_structural_element(*_system, *_assembly, elem, p).release();
}


void MAST::DisplacementOutput::get_weights(std::vector<Real>& w)
{
    for (uint i=0; i<_w.size(); i++)
    {
        w[i] = _w[i];
    }
}
