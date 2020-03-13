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
#include "property_cards/solid_1d_section_element_property_card.h"
#include "elasticity/structural_element_1d.h"
#include "elasticity/structural_system_initialization.h"
#include "base/nonlinear_implicit_assembly.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "base/nonlinear_system.h"
#include "mesh/geom_elem.h"

// Test includes
#include "catch.hpp"
#include "test_helpers.h"
#include "element/structural/1D/mast_structural_element_1d.h"

extern libMesh::LibMeshInit* p_global_init;


TEST_CASE("edge2_linear_extension_structural",
          "[1D][structural][edge][edge2][linear]")
{
    RealMatrixX coords = RealMatrixX::Zero(3, 2);
    coords << -1.0, 1.0, 0.0,
               0.0, 0.0, 0.0;
    TEST::TestStructuralSingleElement1D test_elem(libMesh::EDGE2, coords);

    const Real V0 = test_elem.reference_elem->volume();

    // Set shear coefficient to zero to disable transverse shear stiffness
    test_elem.kappa_zz = 0.0;
    test_elem.kappa_yy = 0.0;

    // Set the offset to zero to disable extension-bending coupling
    test_elem.offset_y = 0.0;
    test_elem.offset_z = 0.0;

    // Set the bending operator to no_bending to disable bending stiffness. This also disables
    // transverse shear stiffness.
    test_elem.section.set_bending_model(MAST::NO_BENDING);

    // Update residual and Jacobian storage since we have modified baseline test element properties.
    test_elem.update_residual_and_jacobian0();

    double val_margin = (test_elem.jacobian0.array().abs()).mean() * 1.490116119384766e-08;

    libMesh::out << "J =\n" << test_elem.jacobian0 << std::endl;

    
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


    SECTION("Internal Jacobian (stiffness matrix) eigenvalue check")
    {
        /*
         * Number of zero eigenvalues should equal the number of rigid body
         * modes.  For 1D extension (including torsion), we have 1 rigid
         * translations along the element's x-axis and 1 rigid rotation about
         * the element's x-axis, for a total of 2 rigid body modes.
         *
         * Note that the use of reduced integration can result in more rigid
         * body modes than expected.
         */
        SelfAdjointEigenSolver<RealMatrixX> eigensolver(test_elem.jacobian0, false);
        RealVectorX eigenvalues = eigensolver.eigenvalues();
        uint nz = 0;
        for (uint i=0; i<eigenvalues.size(); i++)
        {
            if (std::abs(eigenvalues(i))<1e-10)
            {
                nz++;
            }
        }
        REQUIRE(nz == 2);

        /*
         * All non-zero eigenvalues should be positive.
         */
        REQUIRE(eigenvalues.minCoeff()>(-1e-10));
    }


    SECTION("Internal Jacobian (stiffness matrix) is invariant to displacement solution")
    {
        RealVectorX elem_sol = RealVectorX::Zero(test_elem.n_dofs);
        elem_sol << 0.05727841,  0.08896581,  0.09541619, -0.03774913,
                    0.07510557, -0.07122266, -0.00979117, -0.08300009,
                   -0.03453369, -0.05487761, -0.01407677, -0.09268421;
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


    SECTION("Internal Jacobian (stiffness matrix) checks for element aligned along y-axis")
    {
        /*
         * NOTE: We could try to use the transform_element method here, but the
         * issue is that if the sin and cos calculations are not exact, then we
         * may not be perfectly aligned along the y axis like we want.
         */
        RealMatrixX new_coordinates = RealMatrixX::Zero(3,test_elem.n_nodes);
        new_coordinates << 0.0, 0.0, -1.0, 1.0, 0.0, 0.0;
        test_elem.update_coordinates(new_coordinates);
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


    SECTION("Internal Jacobian (stiffness matrix) checks for element aligned along z-axis")
    {
        /*
         * NOTE: We could try to use the transform_element method here, but the
         * issue is that if the sin and cos calculations are not exact, then we
         * may not be perfectly aligned along the z axis like we want.
         */
        RealMatrixX new_coordinates = RealMatrixX::Zero(3,test_elem.n_nodes);
        new_coordinates << 0.0, 0.0, 0.0, 0.0, -1.0, 1.0;
        test_elem.update_coordinates(new_coordinates);
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


    SECTION("Internal Jacobian (stiffness matrix) checks for element arbitrarily scaled, stretched, and rotated")
    {
        // Apply arbitrary transformation to the element
        TEST::transform_element(test_elem.mesh, coords, -5.0, 7.8, -13.1,
                                2.7, 6.4, 20.0, 47.8, -70.1);
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
        TEST::transform_element(test_elem.mesh, coords, 4.1, -6.3, 7.5,
                        4.2, 1.5, -18.0, -24.8, 30.1);

        // Apply arbitrary displacement to the element
        RealVectorX elem_sol = RealVectorX::Zero(test_elem.n_dofs);
        elem_sol << 0.08158724,  0.07991906, -0.00719128,  0.02025461,
                   -0.04602193,  0.05280159,  0.03700081,  0.04636344,
                    0.05559377,  0.06448206,  0.08919238, -0.03079122;
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
