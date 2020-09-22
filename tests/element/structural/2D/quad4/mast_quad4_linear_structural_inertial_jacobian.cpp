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


TEST_CASE("quad4_linear_structural_inertial_consistent", 
          "[quad],[quad4],[dynamic],[2D],[element]")
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
    
    // Set the bending operator
    test_elem.section.set_bending_model(MAST::MINDLIN);

    // Set mass matrix to consistent & compute inertial residual and Jacobians.
    test_elem.section.set_diagonal_mass_matrix(false);
    REQUIRE_FALSE(test_elem.section.if_diagonal_mass_matrix());
    test_elem.update_inertial_residual_and_jacobian0();

    double val_margin = (test_elem.jacobian_xddot0.array().abs()).mean() * 1.490116119384766e-08;

    libMesh::out << "Jac_xddot0:\n" << test_elem.jacobian_xddot0 << std::endl;


    SECTION("Inertial Jacobian (mass matrix) finite difference check")
    {
        TEST::approximate_inertial_jacobian_with_finite_difference(*test_elem.elem,
                                                                   test_elem.elem_solution, test_elem.jacobian_fd);

        val_margin = (test_elem.jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_xddot0),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));
    }


    SECTION("Inertial Jacobian (mass matrix) should be symmetric")
    {
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_xddot0),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_xddot0.transpose())));
    }


    SECTION("Determinant of undeformed inertial Jacobian (mass matrix) should be zero")
    {
        REQUIRE(test_elem.jacobian_xddot0.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("Inertial Jacobian (mass matrix) eigenvalues should all be positive")
    {
        SelfAdjointEigenSolver<RealMatrixX> eigensolver(test_elem.jacobian_xddot0, false);
        RealVectorX eigenvalues = eigensolver.eigenvalues();
        libMesh::out << "Eigenvalues are:\n" << eigenvalues << std::endl;
        REQUIRE(eigenvalues.minCoeff()>0.0);
    }
}


TEST_CASE("quad4_linear_structural_inertial_lumped", 
          "[quad],[quad4],[dynamic],[2D],[element]")
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

    // Set the bending operator
    test_elem.section.set_bending_model(MAST::MINDLIN);

    // Set mass matrix to lumped & compute inertial residual and Jacobians.
    test_elem.section.set_diagonal_mass_matrix(true);
    REQUIRE(test_elem.section.if_diagonal_mass_matrix());
    test_elem.update_inertial_residual_and_jacobian0();

    double val_margin = (test_elem.jacobian_xddot0.array().abs()).mean() * 1.490116119384766e-08;

    libMesh::out << "Jac_xddot0:\n" << test_elem.jacobian_xddot0 << std::endl;


    SECTION("Inertial Jacobian (mass matrix) finite difference check")
    {
        TEST::approximate_inertial_jacobian_with_finite_difference(*test_elem.elem,
                                                                   test_elem.elem_solution, test_elem.jacobian_fd);

        val_margin = (test_elem.jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;

        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_xddot0),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_fd)).margin(val_margin));
    }


    SECTION("Inertial Jacobian (mass matrix) should be symmetric")
    {
        REQUIRE_THAT(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_xddot0),
                     Catch::Approx<double>(TEST::eigen_matrix_to_std_vector(test_elem.jacobian_xddot0.transpose())));
    }


    SECTION("Determinant of undeformed inertial Jacobian (mass matrix) should be zero")
    {
        REQUIRE(test_elem.jacobian_xddot0.determinant() == Approx(0.0).margin(1e-06));
    }


    SECTION("Inertial Jacobian (mass matrix) eigenvalues should all be positive")
    {
        SelfAdjointEigenSolver<RealMatrixX> eigensolver(test_elem.jacobian_xddot0, false);
        RealVectorX eigenvalues = eigensolver.eigenvalues();
        libMesh::out << "Eigenvalues are:\n" << eigenvalues << std::endl;
        REQUIRE(eigenvalues.minCoeff()>0.0);
    }
}
