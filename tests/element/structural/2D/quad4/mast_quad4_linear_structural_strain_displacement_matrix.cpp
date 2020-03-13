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

#define protected public

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
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
#include "mesh/fe_base.h"
#include "numerics/fem_operator_matrix.h"

// Test includes
#include "catch.hpp"
#include "test_helpers.h"
#include "element/structural/2D/mast_structural_element_2d.h"

extern libMesh::LibMeshInit* p_global_init;

/**
 * References
 * ----------
 * https://studiumbook.com/properties-of-shape-function-fea/
 * https://www.ccg.msm.cam.ac.uk/images/FEMOR_Lecture_2.pdf
 */
TEST_CASE("quad4_linear_structural_strain_displacement_matrix", 
          "[quad],[quad4],[linear],][structural],[2D],[element]")
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


    std::unique_ptr<MAST::FEBase> fe(test_elem.geom_elem.init_fe(true, false,
            test_elem.section.extra_quadrature_order(test_elem.geom_elem)));

    const uint n_dofs_per_node = 6;
    MAST::FEMOperatorMatrix
        Bmat_lin,
        Bmat_nl_x,
        Bmat_nl_y,
        Bmat_nl_u,
        Bmat_nl_v,
        Bmat_bend,
        Bmat_vk;
    
    const uint n_phi = (unsigned int)fe->get_phi().size();
    const uint n1 = test_elem.elem->n_direct_strain_components();
    const uint n2 = 6*n_phi;
    const uint n3 = test_elem.elem->n_von_karman_strain_components();
    
    RealMatrixX mat_x = RealMatrixX::Zero(3,2);
    RealMatrixX mat_y = RealMatrixX::Zero(3,2);
    
    RealVectorX strain = RealVectorX::Zero(3);
    
    uint qp = 0;
    
    /**
     * elem->initialize_green_lagrange_strain_operator method populates the
     * Bmat_lin, Bmat_nl_x, Bmat_nl_y, Bmat_nl_u, and Bmat_nl_v matrices.
     */
    Bmat_lin.reinit(n1, test_elem.structural_system.n_vars(), n_phi); // three stress-strain components
    Bmat_nl_x.reinit(2, test_elem.structural_system.n_vars(), n_phi);
    Bmat_nl_y.reinit(2, test_elem.structural_system.n_vars(), n_phi);
    Bmat_nl_u.reinit(2, test_elem.structural_system.n_vars(), n_phi);
    Bmat_nl_v.reinit(2, test_elem.structural_system.n_vars(), n_phi);

    test_elem.elem->initialize_green_lagrange_strain_operator(qp, *fe, test_elem.elem_solution,
        strain, mat_x, mat_y, Bmat_lin, Bmat_nl_x, Bmat_nl_y, Bmat_nl_u,
        Bmat_nl_v);
    
    /**
     * std::unique_ptr<MAST::BendingOperator2D> bend;
     * bend->initialize_bending_strain_operator method populates the Bmat_bend
     * matrix. This part only exists if bending exists in the model.
     */
    Bmat_bend.reinit(n1, test_elem.structural_system.n_vars(), n_phi);
    
    /**
     * elem->initialize_von_karman_strain_operator method populates the Bmat_vk
     * matrix. This part only exists if bending exists in the model AND 
     * nonlinear strains exist in the model.
     */
    Bmat_vk.reinit(n3, test_elem.structural_system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    
    const std::vector<std::vector<libMesh::RealVectorValue>>& dphi = fe->get_dphi();
    
    SECTION("Linear in-plane strain displacement matrix size")
    {
        REQUIRE( Bmat_lin.m() == 3 ); // three strains, e_xx, e_yy, e_xy
        REQUIRE( Bmat_lin.n() == test_elem.n_nodes * n_dofs_per_node);
    }
    
    SECTION("Linear in-plane strain displacement matrix values")
    {
        // TODO: This requires the vector_mult method be working correctly. Is that ok?
        
        // First get a RealMatrixX representation of this matrix
        uint m = Bmat_lin.m();
        uint n = Bmat_lin.n();
        RealMatrixX Bmat_lin_mat = RealMatrixX::Zero(m,n);
        for (uint i=0; i<n; i++)
        {
            RealVectorX Ivec = RealVectorX::Zero(n);
            RealVectorX result = RealVectorX::Zero(m);
            Ivec(i) = 1.0;
            Bmat_lin.vector_mult(result, Ivec);
            Bmat_lin_mat.col(i) = result;
        }
        
        /**
         * Now compare the true values to the expected values
         * 
         * Expected format for Bmat_lin is...
         * [dN1dx, dN2dx, dN3dx, dN4dx,   0,     0,     0,     0,   0, ..., 0;
         *    0,     0,     0,     0,   dN1dy, dN2dy, dN3dy, dN4dy, 0, ..., 0;
         *  dN1dy, dN2dy, dN3dy, dN4dy, dN1dx, dN2dx, dN3dx, dN4dx, 0, ..., 0]
         */
        RealMatrixX Bmat_lin_true = RealMatrixX::Zero(m,n);
        for (uint i=0; i<test_elem.n_nodes; i++)
        {
            Bmat_lin_true(0,i)          = dphi[i][qp](0);
            Bmat_lin_true(2,i+test_elem.n_nodes)  = dphi[i][qp](0);
            
            Bmat_lin_true(1,i+test_elem.n_nodes)  = dphi[i][qp](1);
            Bmat_lin_true(2,i)          = dphi[i][qp](1);
        }
                
        std::vector<double> Bmat_lin_test =    TEST::eigen_matrix_to_std_vector(Bmat_lin_mat);
        std::vector<double> Bmat_lin_required = TEST::eigen_matrix_to_std_vector(Bmat_lin_true);
        
        REQUIRE_THAT(Bmat_lin_test, Catch::Approx<double>(Bmat_lin_required));
    }
}
