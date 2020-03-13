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

// We need access to the protected thermal_residual method to test it
// NOTE: Be careful with this, it could cause unexpected problems
#define protected public

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
#include "base/boundary_condition_base.h"

// Test includes
#include "catch.hpp"
#include "test_helpers.h"
#include "element/structural/1D/mast_structural_element_1d.h"

extern libMesh::LibMeshInit* p_global_init;


TEST_CASE("edge2_linear_structural_thermal_jacobian",
          "[1D],[thermoelastic],[edge],[edge2],[linear],[protected]")
{
    RealMatrixX coords = RealMatrixX::Zero(3, 2);
    coords << -1.0, 1.0, 0.0,
               0.0, 0.0, 0.0;
    TEST::TestStructuralSingleElement1D test_elem(libMesh::EDGE2, coords);

    // Define the Uniform Temperature and Uniform Reference Temperature
    MAST::Parameter temperature("T", 400.0);
    MAST::Parameter ref_temperature("T0", 0.0);
    MAST::ConstantFieldFunction temperature_f("temperature", temperature);
    MAST::ConstantFieldFunction ref_temperature_f("ref_temperature", ref_temperature);

    // Setup the temperature change boundary condition
    MAST::BoundaryConditionBase temperature_load(MAST::TEMPERATURE);
    temperature_load.add(temperature_f);
    temperature_load.add(ref_temperature_f);
    test_elem.discipline.add_volume_load(0, temperature_load);

    const Real V0 = test_elem.reference_elem->volume();

    // Calculate residual and jacobian
    test_elem.elem->thermal_residual(true, test_elem.residual, test_elem.jacobian0, temperature_load);
    
    double val_margin = (test_elem.jacobian0.array().abs()).mean() * 1.490116119384766e-08;
    
    libMesh::out << "R=\n" << test_elem.residual << std::endl;
    libMesh::out << "J =\n" << test_elem.jacobian0 << std::endl;
    
    SECTION("thermal_jacobian_is_zero_matrix")
    {
        // Approximate Jacobian with Finite Difference
        RealMatrixX zero_matrix = RealMatrixX::Zero(test_elem.n_dofs, test_elem.n_dofs);
        
        std::vector<double> test =  TEST::eigen_matrix_to_std_vector(test_elem.jacobian0);
        std::vector<double> truth = TEST::eigen_matrix_to_std_vector(zero_matrix);
        
        REQUIRE_THAT( test, Catch::Approx<double>(truth).margin(0.0) );
    }


    SECTION("Thermoelastic Jacobian finite difference check")
    {
        TEST::approximate_thermal_jacobian_with_finite_difference(*test_elem.elem,
                                                                  test_elem.elem_solution, test_elem.jacobian_fd,
                                                                  temperature_load);

        val_margin = 1.0e-5;

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));
    }


    SECTION("Thermoelastic Jacobian should be symmetric")
    {
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0.transpose())));
    }


    SECTION("Determinant of undeformed thermoelastic Jacobian should be zero")
    {
        REQUIRE(test_elem.jacobian0.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("Thermoelastic Jacobian eigenvalue check")
    {
        /*
         * Linear thermoelastic Jacobian should be independent of the
         * displacements and thus should be a zero matrix.
         */
        SelfAdjointEigenSolver<RealMatrixX> eigensolver(test_elem.jacobian0, false);
        RealVectorX eigenvalues = eigensolver.eigenvalues();
        libMesh::out << "Eigenvalues are:\n" << eigenvalues << std::endl;
        uint nz = 0;
        for (uint i=0; i<eigenvalues.size(); i++)
        {
            if (std::abs(eigenvalues(i))<0.0001220703125)
            {
                nz++;
            }
        }
        REQUIRE( nz == 12);
    }


    SECTION("Thermoelastic Jacobian is invariant to section orientation")
    {
        test_elem.section.clear();
        RealVectorX orientation = RealVectorX::Zero(3);
        orientation(2) = 1.0;
        test_elem.section.y_vector() = orientation;
        test_elem.section.init();
        test_elem.elem->thermal_residual(true, test_elem.residual, test_elem.jacobian, temperature_load);

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0)).margin(val_margin));
    }


    SECTION("Thermoelastic Jacobian is invariant to displacement solution")
    {
        RealVectorX elem_sol = RealVectorX::Zero(test_elem.n_dofs);
        elem_sol << 0.05727841,  0.08896581,  0.09541619, -0.03774913,
                0.07510557, -0.07122266, -0.00979117, -0.08300009,
                -0.03453369, -0.05487761, -0.01407677, -0.09268421;
        test_elem.elem->set_solution(elem_sol);
        test_elem.elem->thermal_residual(true, test_elem.residual, test_elem.jacobian, temperature_load);


        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0)).margin(val_margin));
    }


    SECTION("Thermoelastic Jacobian is invariant to element x-location")
    {
        TEST::transform_element(test_elem.mesh, coords,5.2, 0.0, 0.0,
                                1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.elem->thermal_residual(true, test_elem.residual, test_elem.jacobian, temperature_load);

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0)).margin(val_margin));
    }


    SECTION("Thermoelastic Jacobian is invariant to element y-location")
    {
        TEST::transform_element(test_elem.mesh, coords, 0.0, -11.5, 0.0,
                                1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.elem->thermal_residual(true, test_elem.residual, test_elem.jacobian, temperature_load);

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0)).margin(val_margin));
    }


    SECTION("Thermoelastic Jacobian is invariant to element z-location")
    {
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 7.6,
                                1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.elem->thermal_residual(true, test_elem.residual, test_elem.jacobian, temperature_load);

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian0)).margin(val_margin));
    }


    SECTION("Thermoelastic Jacobian checks for element aligned along y-axis")
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
        test_elem.elem->thermal_residual(true, test_elem.residual, test_elem.jacobian, temperature_load);

        // Finite difference Jacobian check
        TEST::approximate_thermal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd, temperature_load);
        val_margin = 1.0e-5;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }

    SECTION("Thermoelastic Jacobian checks for element aligned along z-axis")
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
        test_elem.elem->thermal_residual(true, test_elem.residual, test_elem.jacobian, temperature_load);

        // Finite difference Jacobian check
        TEST::approximate_thermal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd, temperature_load);
        val_margin = 1.0e-5;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }
    SECTION("Thermoelastic Jacobian checks for element rotated about z-axis")
    {
        // Rotated 63.4 about z-axis at element's centroid
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0,
                                1.0, 1.0, 0.0, 0.0, 63.4);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.elem->thermal_residual(true, test_elem.residual, test_elem.jacobian, temperature_load);

        // Finite difference Jacobian check
        TEST::approximate_thermal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd, temperature_load);
        val_margin = 1.0e-6;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("Thermoelastic Jacobian checks for element rotated about y-axis")
    {
        // Rotated 35.8 about y-axis at element's centroid
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0,
                                1.0, 1.0, 0.0, 35.8, 0.0);
        REQUIRE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.elem->thermal_residual(true, test_elem.residual, test_elem.jacobian, temperature_load);

        // Finite difference Jacobian check
        TEST::approximate_thermal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd, temperature_load);
        val_margin = 1.0e-6;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }
    SECTION("Thermoelastic Jacobian checks for element stretched in x-direction")
    {
        TEST::transform_element(test_elem.mesh, coords, 0.0, 0.0, 0.0,
                                3.2, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        REQUIRE_FALSE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.elem->thermal_residual(true, test_elem.residual, test_elem.jacobian, temperature_load);

        // Finite difference Jacobian check
        TEST::approximate_thermal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd, temperature_load);
        val_margin = 1.0e-6;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("Thermoelastic Jacobian checks for element arbitrarily scaled, stretched, and rotated")
    {
        // Apply arbitrary transformation to the element
        TEST::transform_element(test_elem.mesh, coords, -5.0, 7.8, -13.1,
                                2.7, 6.4, 20.0, 47.8, -70.1);
        REQUIRE_FALSE(test_elem.reference_elem->volume() == Approx(V0));
        test_elem.elem->thermal_residual(true, test_elem.residual, test_elem.jacobian, temperature_load);

        // Finite difference Jacobian check
        TEST::approximate_thermal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd, temperature_load);
        val_margin = 1.0e-6;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("Thermoelastic Jacobian checks for element arbitrarily scaled, stretched, rotated, and displaced")
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
        test_elem.elem->thermal_residual(true, test_elem.residual, test_elem.jacobian, temperature_load);

        // Finite difference Jacobian check
        TEST::approximate_thermal_jacobian_with_finite_difference(*test_elem.elem, test_elem.elem_solution, test_elem.jacobian_fd, temperature_load);
        val_margin = 1.0e-6;
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));

        // Symmetry check
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian.transpose())));

        // Determinant check
        REQUIRE(test_elem.jacobian.determinant() == Approx(0.0).margin(1e-06));
    }
}
