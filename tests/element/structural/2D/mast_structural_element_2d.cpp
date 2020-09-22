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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/elem.h"
#include "libmesh/dof_map.h"

// MAST includes
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "elasticity/structural_element_2d.h"
#include "elasticity/structural_system_initialization.h"
#include "base/nonlinear_implicit_assembly.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "base/nonlinear_system.h"
#include "mesh/geom_elem.h"

// Test includes
#include "catch.hpp"
#include "test_helpers.h"
#include "element/structural/2D/mast_structural_element_2d.h"

extern libMesh::LibMeshInit* p_global_init;

TEST_CASE("structural_element_2d_base_tests",
          "[2D][structural][base]")
{
    RealMatrixX coords = RealMatrixX::Zero(3,4);
    coords << -1.0,  1.0, 1.0, -1.0,
              -1.0, -1.0, 1.0,  1.0,
               0.0,  0.0, 0.0,  0.0;
    TEST::TestStructuralSingleElement2D test_elem(libMesh::QUAD4, coords);

    REQUIRE(test_elem.reference_elem->volume() == TEST::get_shoelace_area(coords));

    SECTION("Element returns proper number of strain components")
    {
        REQUIRE(test_elem.elem->n_direct_strain_components() == 3);
        REQUIRE(test_elem.elem->n_von_karman_strain_components() == 2);
    }
    
    SECTION("Incompatible modes flag returns false for basic element")
    {
        REQUIRE_FALSE(test_elem.elem->if_incompatible_modes());
    }

    SECTION("Isotropic flag returns true for basic element")
    {
        const MAST::ElementPropertyCardBase& elem_section = test_elem.elem->elem_property();
        CHECK( elem_section.if_isotropic() );
    }

    SECTION("Check setting/getting element local solution")
    {
        const libMesh::DofMap& dof_map = test_elem.assembly.system().get_dof_map();
        std::vector<libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices (test_elem.reference_elem, dof_indices);
        uint n_dofs = uint(dof_indices.size());

        RealVectorX elem_solution = 5.3*RealVectorX::Ones(n_dofs);
        test_elem.elem->set_solution(elem_solution);

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(elem_solution),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.elem->local_solution())));
    }

    SECTION("Check setting/getting element local sensitivity solution")
    {
        const libMesh::DofMap& dof_map = test_elem.assembly.system().get_dof_map();;
        std::vector<libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices (test_elem.reference_elem, dof_indices);
        uint n_dofs = uint(dof_indices.size());

        RealVectorX elem_solution_sens = 3.1*RealVectorX::Ones(n_dofs);
        test_elem.elem->set_solution(elem_solution_sens, true);

        const RealVectorX& local_solution_sens = test_elem.elem->local_solution(true);

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(elem_solution_sens),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.elem->local_solution(true))));
    }

    SECTION("Element shape can be transformed")
    {
        const Real V0 = test_elem.reference_elem->volume();

        // Stretch in x-direction
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0, 3.1, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == 12.4);

        // Stretch in y-direction
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0, 1.0, 3.1, 0.0, 0.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == 12.4);

        // Rotation about z-axis
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 60.0);
        REQUIRE(test_elem.reference_elem->volume() == V0);

        // Rotation about y-axis
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 30.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == V0);

        // Rotation about x-axis
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0, 1.0, 1.0, 20.0, 0.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == V0);

        // Shifted in x-direction
        TEST::transform_element(test_elem.mesh, coords, 10.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == V0);

        // Shifted in y-direction
        TEST::transform_element(test_elem.mesh, coords, 0.0, 7.5, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == V0);

        // Shifted in z-direction
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 4.2, 1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == V0);

        // Shear in x
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 4.2, 1.0, 1.0, 0.0, 0.0, 0.0, 5.2, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));

        // Shear in y
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 4.2, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -6.4);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
    }
}
