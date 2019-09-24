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
#include "elasticity/fluid_structure_assembly_elem_operations.h"
#include "elasticity/structural_element_base.h"
#include "property_cards/element_property_card_1D.h"
#include "base/system_initialization.h"
#include "base/physics_discipline_base.h"
#include "mesh/geom_elem.h"


MAST::FluidStructureAssemblyElemOperations::
FluidStructureAssemblyElemOperations():
MAST::AssemblyElemOperations(),
_base_sol       (false),
_qty_type       (MAST::INVALID_QTY) {
    
}



MAST::FluidStructureAssemblyElemOperations::~FluidStructureAssemblyElemOperations() {
    
}


void
MAST::FluidStructureAssemblyElemOperations::
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
MAST::FluidStructureAssemblyElemOperations::init(const MAST::GeomElem& elem) {
    
    libmesh_assert(!_physics_elem);
    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    _physics_elem =
    MAST::build_structural_element(*_system, *_assembly, elem, p).release();
}



void
MAST::FluidStructureAssemblyElemOperations::elem_calculations(bool if_jac,
                                                              RealVectorX& vec,
                                                              RealMatrixX& mat) {
    
    libmesh_assert(if_jac);
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    
    RealMatrixX
    dummy = RealMatrixX::Zero(mat.rows(), mat.cols());
    mat.setZero();
    vec.setZero();
    
    switch (_qty_type) {
        case MAST::MASS: {
            
            e.inertial_residual(true, vec, mat, dummy, dummy);
        }
            break;
            
            
        case MAST::DAMPING: {
            
            e.inertial_residual(true, vec, dummy, mat, dummy);
            e.side_external_residual(true,
                                     vec,
                                     mat,
                                     dummy,
                                     _discipline->side_loads());
            e.volume_external_residual(true,
                                       vec,
                                       mat,
                                       dummy,
                                       _discipline->volume_loads());
        }
            break;
            
        case MAST::STIFFNESS: {
            
            // create
            
            e.internal_residual(true, vec, mat);
            e.inertial_residual(true, vec, dummy, dummy, mat);
            e.side_external_residual(true,
                                     vec,
                                     dummy,
                                     mat,
                                     _discipline->side_loads());
            e.volume_external_residual(true,
                                       vec,
                                       dummy,
                                       mat,
                                       _discipline->volume_loads());
        }
            break;
            
        default:
            libmesh_error(); // should not get here
    }
    
    
}


void
MAST::FluidStructureAssemblyElemOperations::elem_aerodynamic_force_calculations(ComplexVectorX& vec) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    
    ComplexMatrixX
    dummy = ComplexMatrixX::Zero(vec.size(), vec.size());
    vec.setZero();
    
    e.linearized_frequency_domain_side_external_residual(false,
                                                         vec,
                                                         dummy,
                                                         _discipline->side_loads());
    e.linearized_frequency_domain_volume_external_residual(false,
                                                           vec,
                                                           dummy,
                                                           _discipline->volume_loads());
}




void
MAST::FluidStructureAssemblyElemOperations::elem_sensitivity_calculations(const MAST::FunctionBase& f,
                                                                          bool if_jac,
                                                                          RealVectorX& vec,
                                                                          RealMatrixX& mat) {
    
    libmesh_assert(if_jac);
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    
    RealMatrixX
    dummy = RealMatrixX::Zero(vec.size(), vec.size());
    mat.setZero();
    vec.setZero();
    
    switch (_qty_type) {
        case MAST::MASS: {
            
            e.inertial_residual_sensitivity(f, true, vec, mat, dummy, dummy);
        }
            break;
            
        case MAST::DAMPING: {
            
            e.inertial_residual_sensitivity(f, true, vec, dummy, mat, dummy);
            e.side_external_residual_sensitivity(f,
                                                 true,
                                                 vec,
                                                 mat,
                                                 dummy,
                                                 _discipline->side_loads());
            e.volume_external_residual_sensitivity(f,
                                                   true,
                                                   vec,
                                                   mat,
                                                   dummy,
                                                   _discipline->volume_loads());
        }
            break;
            
            
        case MAST::STIFFNESS: {
            
            e.internal_residual_sensitivity(f, true, vec, mat);
            
            // if the linearization is about a base state, then the sensitivity of
            // the base state will influence the sensitivity of the Jacobian
            if (_base_sol)
                e.internal_residual_jac_dot_state_sensitivity(mat);
            
            e.inertial_residual_sensitivity(f, true, vec, dummy, dummy, mat);
            e.side_external_residual_sensitivity(f, true,
                                                 vec,
                                                 dummy,
                                                 mat,
                                                 _discipline->side_loads());
            e.volume_external_residual_sensitivity(f, true,
                                                   vec,
                                                   dummy,
                                                   mat,
                                                   _discipline->volume_loads());
        }
            break;
            
        default:
            libmesh_error(); // should not get here
    }
    
    
}



void
MAST::FluidStructureAssemblyElemOperations::
elem_second_derivative_dot_solution_assembly(RealMatrixX& m) {
    
    libmesh_error(); // to be implemented
}


