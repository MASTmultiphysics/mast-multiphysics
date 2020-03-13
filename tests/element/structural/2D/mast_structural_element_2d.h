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

#ifndef MAST_MAST_STRUCTURAL_ELEMENT_2D_H
#define MAST_MAST_STRUCTURAL_ELEMENT_2D_H

// Test includes
#include "base/mast_mesh.h"

extern libMesh::LibMeshInit* p_global_init;

namespace TEST {

    /**
     *
     */
    class TestStructuralSingleElement2D: public TEST::TestMeshSingleElement {
    public:
        // Material properties.
        MAST::Parameter E;                ///< Modulus of Elasticity
        MAST::Parameter nu;               ///<  Poisson's ratio
        MAST::Parameter rho;              ///<  Density
        MAST::Parameter alpha;            ///<  Coefficient of thermal expansion
        MAST::Parameter cp;               ///<  Specific heat capacity
        MAST::Parameter k;                ///<  Thermal conductivity
        // Section properties.
        MAST::Parameter thickness;  ///< Section thickness
        MAST::Parameter offset;     ///< Section offset
        MAST::Parameter kappa;      ///< Shear coefficient
        // Field functions to distribution parameters throughout model.
        MAST::ConstantFieldFunction E_f;
        MAST::ConstantFieldFunction nu_f;
        MAST::ConstantFieldFunction rho_f;
        MAST::ConstantFieldFunction alpha_f;
        MAST::ConstantFieldFunction cp_f;
        MAST::ConstantFieldFunction k_f;
        MAST::ConstantFieldFunction thickness_f;
        MAST::ConstantFieldFunction offset_f;
        MAST::ConstantFieldFunction kappa_f;
        // Material and property cards.
        MAST::IsotropicMaterialPropertyCard material;
        MAST::Solid2DSectionElementPropertyCard section;
        // Data associated with finite element system/assembly and discipline.
        libMesh::EquationSystems equation_systems;
        MAST::NonlinearSystem& system;
        libMesh::FEType fetype;
        MAST::StructuralSystemInitialization structural_system;
        MAST::PhysicsDisciplineBase discipline;
        MAST::NonlinearImplicitAssembly assembly;
        // Data for actual structural element.
        MAST::GeomElem geom_elem;
        std::unique_ptr<MAST::StructuralElementBase> elem_base;
        MAST::StructuralElement2D* elem;
        // Quick reference to element degree of freedom IDs.
        std::vector<libMesh::dof_id_type> dof_indices;
        uint n_dofs;
        // Element states.
        RealVectorX elem_solution;
        RealVectorX elem_accel;
        // Storage for element residual & Jacobian for zero solution.
        RealVectorX residual;    ///< Vector storage for element's residual vector.
        RealMatrixX jacobian0;   ///< Matrix storage for Jacobian of baseline/undeformed element.
        RealMatrixX jacobian_xdot0;   ///< Matrix storage for velocity Jacobian of baseline/undeformed element.
        RealMatrixX jacobian_xddot0;  ///< Matrix storage for acceleration Jacobian of baseline/undeformed element.
        RealMatrixX jacobian;    ///< Matrix storage for Jacobian of the element in a perturbed/modified state.
        RealMatrixX jacobian_fd; ///< Matrix storage for element Jacobian approximated by finite difference.

        TestStructuralSingleElement2D(libMesh::ElemType e_type, RealMatrixX& coordinates):
        TestMeshSingleElement(e_type, coordinates),
        // Initialize material properties to some default values.
        E("E_param", 72.0e9),
        nu("nu_param", 0.33),
        rho("rho_param", 1420.5),
        alpha("alpha_param", 5.43e-05),
        cp("cp_param",   908.0),
        k("k_param",     237.0),
        // Initialize section properties to some default values.
        thickness("th_param", 0.06),
        offset("off_param", 0.03),
        kappa("kappa_param", 5.0/6.0),
        // Initialize field functions using parameters.
        E_f("E", E),
        nu_f("nu", nu),
        rho_f("rho", rho),
        alpha_f("alpha_expansion", alpha),
        cp_f("cp", cp),
        k_f("k_th", k),
        thickness_f("h", thickness),
        offset_f("off", offset),
        kappa_f("kappa", kappa),
        // Initialize system/discipline/etc.
        equation_systems(mesh),
        system(equation_systems.add_system<MAST::NonlinearSystem>("structural")),
        fetype(libMesh::FIRST, libMesh::LAGRANGE),
        structural_system(system, system.name(), fetype),
        discipline(equation_systems)
        {
            // Configure material card.
            material.add(E_f);
            material.add(nu_f);
            material.add(rho_f);
            material.add(alpha_f);
            material.add(k_f);
            material.add(cp_f);

            // Configure section property card.
            section.add(thickness_f);
            section.add(offset_f);
            section.add(kappa_f);
            section.set_material(material);
            discipline.set_property_for_subdomain(0, section);

            // Setup finite element system and discipline.
            equation_systems.init();
            assembly.set_discipline_and_system(discipline, structural_system);

            // Create the MAST element from the libMesh reference element.
            geom_elem.init(*reference_elem, structural_system);
            elem_base = MAST::build_structural_element(structural_system, geom_elem, section);
            elem = (dynamic_cast<MAST::StructuralElement2D*>(elem_base.get()));

            // Get element DOFs for easy reference in tests.
            const libMesh::DofMap& dof_map = assembly.system().get_dof_map();
            dof_map.dof_indices (reference_elem, dof_indices);
            n_dofs = uint(dof_indices.size());

            // Set element's initial solution and solution sensitivity to zero.
            elem_solution = RealVectorX::Zero(n_dofs);
            elem->set_solution(elem_solution);
            elem->set_solution(elem_solution, true);

            // Set element's initial acceleration and acceleration sensitivity to zero.
            elem_accel = RealVectorX::Zero(n_dofs);
            elem->set_acceleration(elem_accel);
            elem->set_acceleration(elem_accel, true);

            // Calculate residual and jacobian
            residual = RealVectorX::Zero(n_dofs);
            jacobian0 = RealMatrixX::Zero(n_dofs, n_dofs);
            jacobian_xdot0 = RealMatrixX::Zero(n_dofs, n_dofs);
            jacobian_xddot0 = RealMatrixX::Zero(n_dofs, n_dofs);
            jacobian = RealMatrixX::Zero(n_dofs, n_dofs);
            jacobian_fd = RealMatrixX::Zero(n_dofs, n_dofs);
        }
        void update_residual_and_jacobian0() {elem->internal_residual(true, residual, jacobian0);}
        void update_residual_and_jacobian() {elem->internal_residual(true, residual, jacobian);}
        void update_inertial_residual_and_jacobian0() {elem->inertial_residual(true, residual, jacobian_xddot0, jacobian_xdot0, jacobian0);}
    };

}

#endif //MAST_MAST_STRUCTURAL_ELEMENT_2D_H
