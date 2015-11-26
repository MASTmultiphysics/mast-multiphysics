/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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

#ifndef __mast_plate_bending_with_offset_h__
#define __mast_plate_bending_with_offset_h__


// C++ includes
#include <memory>
#include <vector>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/fe_type.h"
#include "libmesh/dof_map.h"



namespace MAST {
    
    // Forward declerations
    class StructuralSystemInitialization;
    class StructuralDiscipline;
    class Parameter;
    class ConstantFieldFunction;
    class IsotropicMaterialPropertyCard;
    class Solid2DSectionElementPropertyCard;
    class DirichletBoundaryCondition;
    class BoundaryConditionBase;
    class StressStrainOutputBase;
    class SectionOffset;
    
    
    struct PlateBendingWithOffset {
        
        
        PlateBendingWithOffset();
        
        
        ~PlateBendingWithOffset();
        
        
        /*!
         *   initializes the object for specified characteristics
         */
        void init(libMesh::ElemType e_type, bool if_vk);

        
        /*!
         *   @returns a pointer to the parameter of the specified name.
         *   If no parameter exists by the specified name, then a \p NULL
         *   pointer is returned and a message is printed with a valid list
         *   of parameters.
         */
        MAST::Parameter* get_parameter(const std::string& nm);
        
        /*!
         *  solves the system and returns the final solution
         */
        const libMesh::NumericVector<Real>&
        solve(bool if_write_output = false);
        
        
        /*!
         *  solves the sensitivity of system and returns the final solution
         */
        const libMesh::NumericVector<Real>&
        sensitivity_solve(MAST::Parameter& p,
                          bool if_write_output = false);
        
        
        /*!
         *   clears the stress data structures for a followup analysis
         */
        void clear_stresss();
        
        bool _initialized;

        // length of domain
        Real _length;

        // width of domain
        Real _width;

        // create the mesh
        libMesh::SerialMesh*           _mesh;
        
        // create the equation system
        libMesh::EquationSystems*      _eq_sys;
        
        // create the libmesh system
        libMesh::NonlinearImplicitSystem*  _sys;
        
        // initialize the system to the right set of variables
        MAST::StructuralSystemInitialization* _structural_sys;
        MAST::StructuralDiscipline*           _discipline;
        
        // create the property functions and add them to the
        MAST::Parameter
        *_th,
        *_E,
        *_nu,
        *_kappa,
        *_press,
        *_zero;
        
        MAST::ConstantFieldFunction
        *_th_f,
        *_E_f,
        *_nu_f,
        *_kappa_f,
        *_press_f;
        
        // Section offset
        MAST::SectionOffset*                       _hoff_f;

        // create the material property card
        MAST::IsotropicMaterialPropertyCard*     _m_card;
        
        // create the element property card
        MAST::Solid2DSectionElementPropertyCard* _p_card;
        
        // create the Dirichlet boundary condition on left edge
        MAST::DirichletBoundaryCondition*     _dirichlet_left;
        
        // create the Dirichlet boundary condition on right edge
        MAST::DirichletBoundaryCondition*     _dirichlet_right;

        // create the Dirichlet boundary condition on bottom edge
        MAST::DirichletBoundaryCondition*     _dirichlet_bottom;

        // create the Dirichlet boundary condition on top edge
        MAST::DirichletBoundaryCondition*     _dirichlet_top;

        
        // create the pressure boundary condition
        MAST::BoundaryConditionBase*             _p_load;
        
        // vector of parameters to evaluate sensitivity wrt
        std::vector<MAST::Parameter*> _params_for_sensitivity;
        
        // output quantity objects to evaluate stress
        std::vector<MAST::StressStrainOutputBase*>  _outputs;
    };
}



#endif //  __mast_plate_bending_with_offset_h__

