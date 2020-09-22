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
#include "elasticity/ks_stress_output.h"
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

// libMesh includes.
#include "libmesh/parallel.h"

MAST::KSStressStrainOutput::KSStressStrainOutput():
MAST::StressStrainOutputBase() {
    
}




MAST::KSStressStrainOutput::~KSStressStrainOutput() {
    
}




void
MAST::KSStressStrainOutput::functional_for_all_elems() {
    
    libmesh_assert(!_if_stress_plot_mode);
    libmesh_assert(!_primal_data_initialized);
    libmesh_assert_greater(_sigma0, 0.);
    
    Real
    e_val            = 0.,
    JxW              = 0.;
    
    _JxW_val         = 0.;
    _sigma_vm_int    = 0.;
    _sigma_vm_p_norm = 0.;
    
    // first find the data with the maximum value, to be used for scaling
    
    // iterate over all element data
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    map_it   =  _stress_data.begin(),
    map_end  =  _stress_data.end();
    
    for ( ; map_it != map_end; map_it++) {
        
        std::vector<MAST::StressStrainOutputBase::Data*>::const_iterator
        vec_it   = map_it->second.begin(),
        vec_end  = map_it->second.end();
        
        for ( ; vec_it != vec_end; vec_it++) {
            
            // ask this data point for the von Mises stress value
            e_val    =   (*vec_it)->von_Mises_stress();
            JxW      =   (*vec_it)->quadrature_point_JxW();
            
            // we do not use absolute value here, since von Mises stress
            // is >= 0.
            _sigma_vm_int  +=  exp(_p_norm_stress * (e_val-_sigma0)/_sigma0) * JxW;
            _JxW_val       +=  JxW;
        }
    }
    
    // sum over all processors, since part of the mesh will exist on the
    // other processors.
    if (!_skip_comm_sum) {
        
        _system->system().comm().sum(_sigma_vm_int);
        _system->system().comm().sum(_JxW_val);
    }
    
    _sigma_vm_p_norm         = 1./_p_norm_stress * log(_sigma_vm_int/_JxW_val);
    _primal_data_initialized = true;
}




void
MAST::KSStressStrainOutput::functional_sensitivity_for_elem
(const MAST::FunctionBase& f,
 const libMesh::dof_id_type e_id,
 Real& dsigma_vm_val_df) const {
    
    libmesh_assert(!_if_stress_plot_mode);
    libmesh_assert(_primal_data_initialized);
    libmesh_assert_greater(_sigma0, 0.);
    
    Real
    e_val         = 0.,
    de_val        = 0.,
    num_sens      = 0.,
    JxW           = 0.;
    
    dsigma_vm_val_df = 0.;
    
    // iterate over all element data
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    map_it   =  _stress_data.find(e_id),
    map_end  =  _stress_data.end();
    
    libmesh_assert(map_it != map_end);
    
    std::vector<MAST::StressStrainOutputBase::Data*>::const_iterator
    vec_it   = map_it->second.begin(),
    vec_end  = map_it->second.end();
    
    for ( ; vec_it != vec_end; vec_it++) {
        
        // ask this data point for the von Mises stress value
        e_val    =   (*vec_it)->von_Mises_stress();
        de_val   =   (*vec_it)->dvon_Mises_stress_dp(f);
        JxW      =   (*vec_it)->quadrature_point_JxW();
        
        num_sens    +=  _p_norm_stress * de_val/_sigma0 * exp(_p_norm_stress * (e_val-_sigma0)/_sigma0) * JxW;
    }
    
    dsigma_vm_val_df = 1./_p_norm_stress / (_sigma_vm_int/_JxW_val) * num_sens / _JxW_val;
}





void
MAST::KSStressStrainOutput::functional_boundary_sensitivity_for_elem
(const MAST::FunctionBase& f,
 const libMesh::dof_id_type e_id,
 Real& dsigma_vm_val_df) const {
    
    libmesh_assert(!_if_stress_plot_mode);
    libmesh_assert(_primal_data_initialized);
    libmesh_assert_greater(_sigma0, 0.);
    
    Real
    e_val         = 0.,
    JxW_Vn        = 0.,
    num_sens      = 0.,
    denom_sens    = 0.;
    
    dsigma_vm_val_df = 0.;
    
    // iterate over all element data
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    map_it   =  _boundary_stress_data.find(e_id),
    map_end  =  _boundary_stress_data.end();
    
    libmesh_assert(map_it != map_end);
    
    std::vector<MAST::StressStrainOutputBase::Data*>::const_iterator
    vec_it   = map_it->second.begin(),
    vec_end  = map_it->second.end();
    
    for ( ; vec_it != vec_end; vec_it++) {
        
        // ask this data point for the von Mises stress value
        e_val    =   (*vec_it)->von_Mises_stress();
        JxW_Vn   =   (*vec_it)->quadrature_point_JxW();
        
        denom_sens  +=  JxW_Vn;
        num_sens    +=  exp(_p_norm_stress * (e_val-_sigma0)/_sigma0) * JxW_Vn;
    }
    
    dsigma_vm_val_df = 1./_p_norm_stress / (_sigma_vm_int/_JxW_val) *
    (num_sens / _JxW_val - _sigma_vm_int / pow(_JxW_val, 2.) * denom_sens);
}




void
MAST::KSStressStrainOutput::functional_state_derivartive_for_elem(const libMesh::dof_id_type e_id,
                                                                  RealVectorX& dq_dX) const {
    libmesh_assert(!_if_stress_plot_mode);
    libmesh_assert(_primal_data_initialized);
    libmesh_assert_greater(_sigma0, 0.);
    
    Real
    e_val         = 0.,
    JxW           = 0.;
    
    RealVectorX
    num_sens       = RealVectorX::Zero(dq_dX.size()),
    de_val         = RealVectorX::Zero(dq_dX.size());
    dq_dX.setZero();
    
    // first find the data with the maximum value, to be used for scaling
    
    // iterate over all element data
    std::map<const libMesh::dof_id_type, std::vector<MAST::StressStrainOutputBase::Data*> >::const_iterator
    map_it   =  _stress_data.find(e_id),
    map_end  =  _stress_data.end();
    
    // make sure that the data exists
    libmesh_assert(map_it != map_end);
    
    std::vector<MAST::StressStrainOutputBase::Data*>::const_iterator
    vec_it   = map_it->second.begin(),
    vec_end  = map_it->second.end();
    
    
    for ( ; vec_it != vec_end; vec_it++) {
        
        // ask this data point for the von Mises stress value
        e_val    =   (*vec_it)->von_Mises_stress();
        de_val   =   (*vec_it)->dvon_Mises_stress_dX();
        JxW      =   (*vec_it)->quadrature_point_JxW();
        
        num_sens    += _p_norm_stress * de_val/_sigma0 * exp(_p_norm_stress * (e_val-_sigma0)/_sigma0) * JxW;
    }
    
    dq_dX = 1./_p_norm_stress / (_sigma_vm_int/_JxW_val) * num_sens / _JxW_val;
}




