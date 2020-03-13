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


TEST_CASE("quad4_linear_structural", 
          "[quad],[quad4],[linear],[structural],[2D],[element]")
{
    RealMatrixX coords = RealMatrixX::Zero(3,4);
    coords << -1.0,  1.0, 1.0, -1.0,
            -1.0, -1.0, 1.0,  1.0,
            0.0,  0.0, 0.0,  0.0;
    TEST::TestStructuralSingleElement2D test_elem(libMesh::QUAD4, coords);

    const Real V0 = test_elem.reference_elem->volume();
    REQUIRE(test_elem.reference_elem->volume() == TEST::get_shoelace_area(coords));

    // Set the strain type to linear for the section
    test_elem.section.set_strain(MAST::LINEAR_STRAIN);

    // Update residual and Jacobian storage since we have modified baseline test element properties.
    test_elem.update_residual_and_jacobian0();

    double val_margin = (test_elem.jacobian0.array().abs()).mean() * 1.490116119384766e-08;


    SECTION("Internal Jacobian (stiffness matrix) finite difference check")
    {
        TEST::approximate_internal_jacobian_with_finite_difference(*test_elem.elem,
                                                                   test_elem.elem_solution, test_elem.jacobian_fd);

        val_margin = (test_elem.jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));
    }


    SECTION("Internal Jacobian (stiffness matrix) should be symmetric")
    {
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0.transpose())));
    }


    SECTION("Determinant of undeformed internal Jacobian (stiffness matrix) should be zero")
    {
        REQUIRE(test_elem.jacobian0.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("Internal Jacobian (stiffness matrix) is invariant to displacement solution")
    {
        RealVectorX elem_sol = RealVectorX::Zero(test_elem.n_dofs);
        elem_sol << -0.04384355,  0.03969142, -0.09470648, -0.05011107,
                -0.02989082, -0.01205296,  0.08846868,  0.04522207,
                0.06435953, -0.07282706,  0.09307561, -0.06250143,
                0.03332844, -0.00040089, -0.00423108, -0.07258241,
                0.06636534, -0.08421098, -0.0705489 , -0.06004976,
                0.03873095, -0.09194373,  0.00055061,  0.046831;
        test_elem.elem->set_solution(elem_sol);
        test_elem.update_residual_and_jacobian();

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0)).margin(val_margin));
    }


    SECTION("Internal Jacobian (stiffness matrix) is invariant to element x-location")
    {
        TEST::transform_element(test_elem.mesh, coords,5.2, 0.0, 0.0,
                                1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.update_residual_and_jacobian();

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0)).margin(val_margin));
    }


    SECTION("Internal Jacobian (stiffness matrix) is invariant to element y-location")
    {
        TEST::transform_element(test_elem.mesh, coords, 0.0, -11.5, 0.0,
                                1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.update_residual_and_jacobian();

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0)).margin(val_margin));
    }


    SECTION("Internal Jacobian (stiffness matrix) is invariant to element z-location")
    {
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 7.6,
                                1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.update_residual_and_jacobian();

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0)).margin(val_margin));
    }


    SECTION("Internal Jacobian (stiffness matrix) checks for element rotated about z-axis")
    {
        // Rotated 63.4 about z-axis at element's centroid
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0,
                                1.0, 1.0, 0.0, 0.0, 63.4);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.update_residual_and_jacobian();

        // Finite difference Jacobian check
        TEST::approximate_internal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd);
        val_margin = (test_elem.jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("Internal Jacobian (stiffness matrix) checks for element rotated about y-axis")
    {
        // Rotated 35.8 about y-axis at element's centroid
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0,
                                1.0, 1.0, 0.0, 35.8, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.update_residual_and_jacobian();

        // Finite difference Jacobian check
        TEST::approximate_internal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd);
        val_margin = (test_elem.jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("Internal Jacobian (stiffness matrix) checks for element rotated about x-axis")
    {
        // Rotated 15.8 about x-axis at element's centroid
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0,
                                1.0, 1.0, 15.8, 0.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.update_residual_and_jacobian();

        // Finite difference Jacobian check
        TEST::approximate_internal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd);
        val_margin = (test_elem.jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("\"Internal Jacobian (stiffness matrix) checks for element sheared in x-direction\"")
    {
        // Shear element in x-direction.
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0,
                                1.0, 1.0, 0.0, 0.0, 0.0, 6.7, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.update_residual_and_jacobian();

        // Finite difference Jacobian check
        TEST::approximate_internal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd);
        val_margin = (test_elem.jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("\"Internal Jacobian (stiffness matrix) checks for element sheared in y-direction\"")
    {
        // Shear element in x-direction.
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0,
                                1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -11.2);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.update_residual_and_jacobian();

        // Finite difference Jacobian check
        TEST::approximate_internal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd);
        val_margin = (test_elem.jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("Internal Jacobian (stiffness matrix) checks for element stretched in x-direction")
    {
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0,
                                3.2, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        REQUIRE_FALSE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.update_residual_and_jacobian();

        // Finite difference Jacobian check
        TEST::approximate_internal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd);
        val_margin = (test_elem.jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }

    SECTION("Internal Jacobian (stiffness matrix) checks for element stretched in y-direction")
    {
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0,
                                1.0, 0.64, 0.0, 0.0, 0.0, 0.0, 0.0);
        REQUIRE_FALSE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.update_residual_and_jacobian();

        // Finite difference Jacobian check
        TEST::approximate_internal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd);
        val_margin = (test_elem.jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("Internal Jacobian (stiffness matrix) checks for element arbitrarily scaled, stretched, rotated, and sheared")
    {
        // Apply arbitrary transformation to the element
        TEST::transform_element(test_elem.mesh, coords, -5.0, 7.8, -13.1,
                                2.7, 6.4, 20.0, 47.8, -70.1, 5.7, -6.3);
        REQUIRE_FALSE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.update_residual_and_jacobian();

        // Finite difference Jacobian check
        TEST::approximate_internal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd);
        val_margin = (test_elem.jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("Internal Jacobian (stiffness matrix) checks for element arbitrarily scaled, stretched, rotated, and displaced")
    {
        // Apply arbitrary transformation to the element
        TEST::transform_element(test_elem.mesh, coords, 4.1, -6.3, 7.5, 4.2, 1.5, -18.0, -24.8,
                                30.1, -3.2, 5.4);

        // Apply arbitrary displacement to the element
        RealVectorX elem_sol = RealVectorX::Zero(test_elem.n_dofs);
        elem_sol << -0.04384355,  0.03969142, -0.09470648, -0.05011107,
                -0.02989082, -0.01205296,  0.08846868,  0.04522207,
                0.06435953, -0.07282706,  0.09307561, -0.06250143,
                0.03332844, -0.00040089, -0.00423108, -0.07258241,
                0.06636534, -0.08421098, -0.0705489 , -0.06004976,
                0.03873095, -0.09194373,  0.00055061,  0.046831;
        test_elem.elem->set_solution(elem_sol);

        REQUIRE_FALSE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.update_residual_and_jacobian();

        // Finite difference Jacobian check
        TEST::approximate_internal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd);
        val_margin = (test_elem.jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }
}
