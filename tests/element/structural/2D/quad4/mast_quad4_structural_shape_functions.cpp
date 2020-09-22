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

// NOTE: Be careful with this, it could cause issues.  Needed to access
// protected members to modify them for finite difference sensitivity check.
#define protected public

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/point.h"
#include "libmesh/elem.h"
#include "libmesh/face_quad4.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"

// MAST includes
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "elasticity/structural_element_2d.h"
#include "elasticity/structural_system_initialization.h"
#include "base/physics_discipline_base.h"
#include "base/nonlinear_implicit_assembly.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "base/nonlinear_system.h"
#include "elasticity/structural_element_base.h"
#include "mesh/geom_elem.h"
#include "mesh/fe_base.h"

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
TEST_CASE("quad4_structural_shape_functions", 
          "[quad],[quad4],[structural],[2D],[element]")
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
    
    SECTION("Number of shape functions equal number of nodes")
    {
        REQUIRE(fe->get_phi().size() == test_elem.n_nodes);
    }


    SECTION("Shape functions sum to one at quadrature points")
    {
        // Get the shape function values at the quadrature points
        const std::vector<std::vector<Real>>& phi = fe->get_phi();
        
        uint n_qps = phi[0].size();
        
        for (uint i=0; i<n_qps; i++) // Iterate Over Quadrature Points
        {
            Real phi_j_sum = 0.0;
            for (uint j=0; j<phi.size(); j++) // Iterative Over Shape Functions
            {
                //libMesh::out << "phi[" << j << "][" << i << "] = " << phi[j][i] << std::endl;
                phi_j_sum += phi[j][i];
            }
            REQUIRE(phi_j_sum == Approx(1.0));
        }
    }


    SECTION("Shape function value is 1 at it's node and zero at other nodes")
    {
        // Redfine the points where the shape functions are calculated
        // Default value is the quadrature points
        std::vector<libMesh::Point> pts;
        pts.reserve(test_elem.n_nodes);
        for (uint i=0; i<test_elem.n_nodes; i++)
        {
            pts.push_back(libMesh::Point(coords(0,i), coords(1,i), coords(2,i)));
        }
        delete fe->_fe;
        fe->_fe = nullptr;
        fe->_initialized = false;
        fe->init(test_elem.geom_elem, true, &pts);

        // Get the shape function values at the node points defined above
        const std::vector<std::vector<Real>>& phi = fe->get_phi();

        uint n_qps = phi[0].size();

        for (uint i=0; i<n_qps; i++) // Iterate Over Node Points
        {
            for (uint j=0; j<phi.size(); j++) // Iterative Over Shape Functions
            {
                //libMesh::out << "phi[" << j << "][" << i << "] = " << phi[j][i] << std::endl;
                if (i==j)
                {
                    REQUIRE(phi[j][i] == Approx(1.0));
                }
                else
                {
                    REQUIRE(phi[j][i] == Approx(0.0));
                }
            }
        }
    }


    SECTION("Shape function derivatives sum to zero at quadrature points")
    {
        // Get the shape function derivative values at the quadrature points
        const std::vector<std::vector<libMesh::RealVectorValue>>& dphi = fe->get_dphi();

        uint n_qps = dphi[0].size();

        for (uint i=0; i<n_qps; i++) // Iterate Over Quadrature Points
        {
            Real dphi_i_dx_sum = 0.0;
            Real dphi_i_dy_sum = 0.0;
            for (uint j=0; j<dphi.size(); j++) // Iterative Over Shape Functions
            {
                dphi_i_dx_sum += dphi[j][i](0);
                dphi_i_dy_sum += dphi[j][i](1);
            }
            REQUIRE(dphi_i_dx_sum == Approx(0.0));
            REQUIRE(dphi_i_dy_sum == Approx(0.0));
        }
    }


    SECTION("Shape function xi derivative finite difference check")
    {
        // Get the shape function derivative values at the quadrature points
        const std::vector<std::vector<libMesh::RealVectorValue>>& dphi = fe->get_dphi();

        uint n_qps = dphi[0].size();

        RealMatrixX dphi_dxi_0 = RealMatrixX::Zero(test_elem.n_nodes, n_qps);
        RealMatrixX phi_xi_h = RealMatrixX::Zero(test_elem.n_nodes, n_qps);
        RealMatrixX phi_xi_n = RealMatrixX::Zero(test_elem.n_nodes, n_qps);

        for (uint i=0; i<n_qps; i++)
        {
            for (uint j=0; j<test_elem.n_nodes; j++)
            {
                dphi_dxi_0(j,i) = dphi[j][i](0);
            }
        }

        // Get the quadrature points
        const std::vector<libMesh::Point>& q_pts = fe->get_qpoints();

        // Shift the Quadrature Points in the xi direction
        Real delta = 0.0001220703125; // sqrt(sqrt(eps))
        std::vector<libMesh::Point> pts;
        pts.reserve(test_elem.n_nodes);

        // Shift quadrature points in +xi direction
        for (uint i=0; i<test_elem.n_nodes; i++)
        {
            pts.push_back(libMesh::Point(q_pts[i](0)+delta, q_pts[i](1), q_pts[i](2)));
        }
        delete fe->_fe;
        fe->_fe = nullptr;
        fe->_initialized = false;
        fe->init(test_elem.geom_elem, true, &pts);
        const std::vector<std::vector<Real>>& phi_xih = fe->get_phi();
        for (uint i=0; i<n_qps; i++)
        {
            for (uint j=0; j<test_elem.n_nodes; j++)
            {
                phi_xi_h(j,i) = phi_xih[j][i];
            }
        }

        // Shift quadrature points in -xi direction
        for (uint i=0; i<test_elem.n_nodes; i++)
        {
            pts[i] = libMesh::Point(q_pts[i](0)-delta, q_pts[i](1), q_pts[i](2));
        }
        delete fe->_fe;
        fe->_fe = nullptr;
        fe->_initialized = false;
        fe->init(test_elem.geom_elem, true, &pts);
        const std::vector<std::vector<Real>>& phi_xin = fe->get_phi();
        for (uint i=0; i<n_qps; i++)
        {
            for (uint j=0; j<test_elem.n_nodes; j++)
            {
                phi_xi_n(j,i) = phi_xin[j][i];
            }
        }

        // Calculate second order central difference approximation to dphi_dxi
        RealMatrixX dphi_dxi_fd = RealMatrixX::Zero(test_elem.n_nodes, n_qps);
        for (uint i=0; i<n_qps; i++) // Iterate Over Quadrature Points
        {
            for (uint j=0; j<test_elem.n_nodes; j++) // Iterative Over Shape Functions
            {
                dphi_dxi_fd(j,i) = (phi_xi_h(j,i) - phi_xi_n(j,i))/(2.0*delta) ;
            }
        }
        //libMesh::out << "dphi_dxi:\n" << dphi_dxi_0 << std::endl;
        //libMesh::out << "dphi_dxi_fd:\n" << dphi_dxi_fd << std::endl;

        std::vector<double> dPhi_dxi =    TEST::eigen_matrix_to_std_vector(dphi_dxi_0);
        std::vector<double> dPhi_dxi_fd = TEST::eigen_matrix_to_std_vector(dphi_dxi_fd);

        REQUIRE_THAT( dPhi_dxi, Catch::Approx<double>(dPhi_dxi_fd) );
    }


    SECTION("Shape function eta derivative finite difference check")
    {
        // Get the shape function derivative values at the quadrature points
        const std::vector<std::vector<libMesh::RealVectorValue>>& dphi = fe->get_dphi();

        uint n_qps = dphi[0].size();

        RealMatrixX dphi_deta_0 = RealMatrixX::Zero(test_elem.n_nodes, n_qps);
        RealMatrixX phi_eta_h = RealMatrixX::Zero(test_elem.n_nodes, n_qps);
        RealMatrixX phi_eta_n = RealMatrixX::Zero(test_elem.n_nodes, n_qps);

        for (uint i=0; i<n_qps; i++)
        {
            for (uint j=0; j<test_elem.n_nodes; j++)
            {
                dphi_deta_0(j,i) = dphi[j][i](1);
            }
        }

        // Get the quadrature points
        const std::vector<libMesh::Point>& q_pts = fe->get_qpoints();

        // Shift the Quadrature Points in the xi direction
        Real delta = 0.0001220703125; // sqrt(sqrt(eps))
        std::vector<libMesh::Point> pts;
        pts.reserve(test_elem.n_nodes);

        // Shift quadrature points in +eta direction
        for (uint i=0; i<test_elem.n_nodes; i++)
        {
            pts.push_back(libMesh::Point(q_pts[i](0), q_pts[i](1)+delta, q_pts[i](2)));
        }
        delete fe->_fe;
        fe->_fe = nullptr;
        fe->_initialized = false;
        fe->init(test_elem.geom_elem, true, &pts);
        const std::vector<std::vector<Real>>& phi_etah = fe->get_phi();
        for (uint i=0; i<n_qps; i++)
        {
            for (uint j=0; j<test_elem.n_nodes; j++)
            {
                phi_eta_h(j,i) = phi_etah[j][i];
            }
        }

        // Shift quadrature points in -eta direction
        for (uint i=0; i<test_elem.n_nodes; i++)
        {
            pts[i] = libMesh::Point(q_pts[i](0), q_pts[i](1)-delta, q_pts[i](2));
        }
        delete fe->_fe;
        fe->_fe = nullptr;
        fe->_initialized = false;
        fe->init(test_elem.geom_elem, true, &pts);
        const std::vector<std::vector<Real>>& phi_etan = fe->get_phi();
        for (uint i=0; i<n_qps; i++)
        {
            for (uint j=0; j<test_elem.n_nodes; j++)
            {
                phi_eta_n(j,i) = phi_etan[j][i];
            }
        }

        // Calculate second order central difference approximation to dphi_dxi
        RealMatrixX dphi_deta_fd = RealMatrixX::Zero(test_elem.n_nodes, n_qps);
        for (uint i=0; i<n_qps; i++) // Iterate Over Quadrature Points
        {
            for (uint j=0; j<test_elem.n_nodes; j++) // Iterative Over Shape Functions
            {
                dphi_deta_fd(j,i) = (phi_eta_h(j,i) - phi_eta_n(j,i))/(2.0*delta) ;
            }
        }
        //libMesh::out << "dphi_deta:\n" << dphi_deta_0 << std::endl;
        //libMesh::out << "dphi_deta_fd:\n" << dphi_deta_fd << std::endl;

        std::vector<double> dPhi_deta =    TEST::eigen_matrix_to_std_vector(dphi_deta_0);
        std::vector<double> dPhi_deta_fd = TEST::eigen_matrix_to_std_vector(dphi_deta_fd);

        REQUIRE_THAT( dPhi_deta, Catch::Approx<double>(dPhi_deta_fd) );
    }
}
