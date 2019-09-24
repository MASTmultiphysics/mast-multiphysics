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
#include "elasticity/stress_temperature_adjoint.h"
#include "elasticity/structural_element_base.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "base/assembly_base.h"
#include "mesh/fe_base.h"
#include "mesh/geom_elem.h"

// libMesh includes
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"

MAST::StressTemperatureAdjoint::StressTemperatureAdjoint(MAST::StressStrainOutputBase& stress):
MAST::StressStrainOutputBase(),
_stress             (stress),
_thermal_assembly   (nullptr) {
    
}



MAST::StressTemperatureAdjoint::~StressTemperatureAdjoint() {
    
}


void
MAST::StressTemperatureAdjoint::
set_thermal_assembly(MAST::AssemblyBase& thermal_assembly) {
    
    libmesh_assert(!_thermal_assembly);
    _thermal_assembly = &thermal_assembly;
}



void
MAST::StressTemperatureAdjoint::
set_structural_solutions(const libMesh::NumericVector<Real>& sol,
                         const libMesh::NumericVector<Real>& adj_sol) {
    
    _structural_sol.reset(_assembly->build_localized_vector(_system->system(),
                                                            sol).release());
    _structural_adjoint.reset(_assembly->build_localized_vector(_system->system(),
                                                                adj_sol).release());
}


void
MAST::StressTemperatureAdjoint::set_elem_solution(const RealVectorX& sol) {
    
    // This solves for the adjoint of stress with dependence on temperature.
    // The sol vector here is the temperature solution since this is called
    // from within the thermal system assembly, but the stress object
    // requires the structural solution. Hence, the sol vector is ignored.
    // Instead, the structural solution is extracted and provided to the
    // element
    
    
    libmesh_assert(_physics_elem);

    MAST::StructuralElementBase* p_elem =
    dynamic_cast<MAST::StructuralElementBase*>(_physics_elem);
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    dof_map.dof_indices (&p_elem->elem().get_reference_elem(), dof_indices);
    
    RealVectorX
    str_sol   = RealVectorX::Zero(dof_indices.size());
    
    // get the element structural solution
    for (unsigned int i=0; i<dof_indices.size(); i++)
        str_sol(i) = (*_structural_sol)(dof_indices[i]);

    _physics_elem->set_solution(str_sol);
}


void
MAST::StressTemperatureAdjoint::output_derivative_for_elem(RealVectorX& dq_dX) {
    
    libmesh_assert(_physics_elem);
    libmesh_assert(!_stress.stress_plot_mode());
    libmesh_assert(_stress.primal_data_initialized());


    
    dq_dX.setZero();
    
    const libMesh::Elem& e = _physics_elem->elem().get_reference_elem();
    
    if (_stress.if_evaluate_for_element(_physics_elem->elem())) {

        MAST::StructuralElementBase* p_elem =
        dynamic_cast<MAST::StructuralElementBase*>(_physics_elem);

        std::vector<libMesh::dof_id_type> dof_indices;
        const libMesh::DofMap& dof_map = _system->system().get_dof_map();
        
        dof_map.dof_indices (&e, dof_indices);
        
        RealVectorX
        str_adj   = RealVectorX::Zero(dof_indices.size());
        
        RealMatrixX
        mat       = RealMatrixX::Zero(dof_indices.size(), dq_dX.size());
        
        // get the element structural adjoint solution
        for (unsigned int i=0; i<dof_indices.size(); i++)
            str_adj(i) = (*_structural_adjoint)(dof_indices[i]);
        
        
        if (_stress.get_thermal_load_for_elem(_physics_elem->elem())) {
            
            std::unique_ptr<MAST::FEBase> fe(_physics_elem->elem().init_fe(true, false));

            p_elem->calculate_stress_temperature_derivative(*fe, _stress);
            p_elem->thermal_residual_temperature_derivative(*fe, mat);
        }
        
        _stress.functional_state_derivartive_for_elem(e.id(), dq_dX);
        dq_dX += str_adj.transpose() * mat;
    }
}



