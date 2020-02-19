// NOTE: Be careful with this, it could cause issues.  Needed to access
// protected members to modify them for finite difference sensitivity check.
#define protected public

#include "catch.hpp"

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


// Custom includes
#include "test_helpers.h"

#define pi 3.14159265358979323846

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
    const int n_elems = 1;
    const int n_nodes = 4;
    
    // Point Coordinates
    RealMatrixX temp = RealMatrixX::Zero(3,4);
    temp << -1.0,  1.0, 1.0, -1.0, 
            -1.0, -1.0, 1.0,  1.0, 
             0.0,  0.0, 0.0,  0.0;
    const RealMatrixX X = temp;
    
    /**
     *  First create the mesh with the one element we are testing.
     */
    // Setup the mesh properties
    libMesh::ReplicatedMesh mesh(p_global_init->comm());
    mesh.set_mesh_dimension(2);
    mesh.set_spatial_dimension(2);
    mesh.reserve_elem(n_elems);
    mesh.reserve_nodes(n_nodes);
    
    // Add nodes to the mesh
    for (uint i=0; i<n_nodes; i++)
    {
        mesh.add_point(libMesh::Point(X(0,i), X(1,i), X(2,i)));
    }
    
    // Add the element to the mesh
    libMesh::Elem *reference_elem = new libMesh::Quad4;
    reference_elem->set_id(0);    
    reference_elem->subdomain_id() = 0;
    reference_elem = mesh.add_elem(reference_elem);
    for (int i=0; i<n_nodes; i++)
    {
        reference_elem->set_node(i) = mesh.node_ptr(i);
    }
    
    // Prepare the mesh for use
    mesh.prepare_for_use();
    //mesh.print_info();
    
    const Real elem_volume = reference_elem->volume();
    // Calculate true volume using 2D shoelace formula
    Real true_volume = get_shoelace_area(X);
    
    // Ensure the libMesh element has the expected volume
    REQUIRE( elem_volume == true_volume );
            
    /**
     *  Setup the material and section properties to be used in the element
     */
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    MAST::Parameter kappa("kappa_param", 5.0/6.0);    // Shear coefficient
    
    // Define Section Properties as MAST Parameters
    MAST::Parameter thickness("th_param", 0.06);      // Section thickness
    MAST::Parameter offset("off_param", 0.03);        // Section offset
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction kappa_f("kappa", kappa);
    MAST::ConstantFieldFunction thickness_f("h", thickness);
    MAST::ConstantFieldFunction offset_f("off", offset);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);                                             
    material.add(nu_f);
    
    // Initialize the section
    MAST::Solid2DSectionElementPropertyCard section;
    
    // Add the section property constant field functions to the section card
    section.add(thickness_f);
    section.add(offset_f);
    section.add(kappa_f);
    
    // Add the material card to the section card
    section.set_material(material);
    
    
    // Set the strain type to linear for the section
    section.set_strain(MAST::LINEAR_STRAIN);
    
    /**
     *  Now we setup the structural system we will be solving.
     */
    libMesh::EquationSystems equation_systems(mesh);
    
    MAST::NonlinearSystem& system = equation_systems.add_system<MAST::NonlinearSystem>("structural");
    
    libMesh::FEType fetype(libMesh::FIRST, libMesh::LAGRANGE);
    
    MAST::StructuralSystemInitialization structural_system(system, 
                                                           system.name(), 
                                                           fetype);
    
    MAST::PhysicsDisciplineBase discipline(equation_systems);
    
    discipline.set_property_for_subdomain(0, section);
    
    equation_systems.init();
    //equation_systems.print_info();
    
    MAST::NonlinearImplicitAssembly assembly;
    assembly.set_discipline_and_system(discipline, structural_system);
    
    // Create the MAST element from the libMesh reference element
    MAST::GeomElem geom_elem;
    geom_elem.init(*reference_elem, structural_system);
    std::unique_ptr<MAST::StructuralElementBase> elem_base = build_structural_element(structural_system, geom_elem, section);
    
    // Cast the base structural element as a 2D structural element
    MAST::StructuralElement2D* elem = (dynamic_cast<MAST::StructuralElement2D*>(elem_base.get()));
    
    // Get element DOFs
    const libMesh::DofMap& dof_map = assembly.system().get_dof_map();
    std::vector<libMesh::dof_id_type> dof_indices;
    dof_map.dof_indices (reference_elem, dof_indices);
    uint n_dofs = uint(dof_indices.size());
    
    // Set element's initial solution and solution sensitivity to zero
    RealVectorX elem_solution = RealVectorX::Zero(n_dofs);
    elem->set_solution(elem_solution);
    elem->set_solution(elem_solution, true);
    
    std::unique_ptr<MAST::FEBase> fe(geom_elem.init_fe(true, false, 
                                     section.extra_quadrature_order(geom_elem)));
    
    
    SECTION("number of shape functions equal number of nodes")
    {
        // Get the shape function values at the quadrature points
        const std::vector<std::vector<Real>>& phi = fe->get_phi();
        REQUIRE( phi.size() == n_nodes );
    }
    
    
    SECTION("shape function summation to one")                   
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
            REQUIRE( phi_j_sum == Approx(1.0) );
        }
    }
    
    
    /**
     * Shape function value at each node is 1 for the shape function 
     * associated with that node and zero for all other nodes.
     */
    SECTION("shape function value is 1 at shape function's node, zero at other nodes")
    {
        // Redfine the points where the shape functions are calculated
        // Default value is the quadrature points
        std::vector<libMesh::Point> pts;
        pts.reserve(n_nodes);
        for (uint i=0; i<n_nodes; i++)
        {
            pts.push_back(libMesh::Point(X(0,i), X(1,i), X(2,i)));
        }
        delete fe->_fe;
        fe->_fe = nullptr;
        fe->_initialized = false;
        fe->init(geom_elem, true, &pts);
        
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
                    REQUIRE( phi[j][i] == Approx(1.0) );
                }
                else
                {
                    REQUIRE( phi[j][i] == Approx(0.0) );
                }
            }
        }
    }
    
    
    SECTION("shape function derivative summation to zero")
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
            REQUIRE( dphi_i_dx_sum == Approx(0.0) );
            REQUIRE( dphi_i_dy_sum == Approx(0.0) );
        }
    }
    

    SECTION("shape function xi derivative finite difference check")
    {
        // Get the shape function derivative values at the quadrature points
        const std::vector<std::vector<libMesh::RealVectorValue>>& dphi = fe->get_dphi();
        
        uint n_qps = dphi[0].size();
        
        RealMatrixX dphi_dxi_0 = RealMatrixX::Zero(n_nodes, n_qps);
        RealMatrixX phi_xi_h = RealMatrixX::Zero(n_nodes, n_qps);
        RealMatrixX phi_xi_n = RealMatrixX::Zero(n_nodes, n_qps);
        
        for (uint i=0; i<n_qps; i++)
        {
            for (uint j=0; j<n_nodes; j++)
            {
                dphi_dxi_0(j,i) = dphi[j][i](0);
            }
        }
        
        // Get the quadrature points
        const std::vector<libMesh::Point>& q_pts = fe->get_qpoints();
        
        // Shift the Quadrature Points in the xi direction
        Real delta = 0.0001220703125; // sqrt(sqrt(eps))
        std::vector<libMesh::Point> pts;
        pts.reserve(n_nodes);
        
        // Shift quadrature points in +xi direction
        for (uint i=0; i<n_nodes; i++)
        {
            pts.push_back(libMesh::Point(q_pts[i](0)+delta, q_pts[i](1), q_pts[i](2)));
        }
        delete fe->_fe;
        fe->_fe = nullptr;
        fe->_initialized = false;
        fe->init(geom_elem, true, &pts);
        const std::vector<std::vector<Real>>& phi_xih = fe->get_phi();
        for (uint i=0; i<n_qps; i++)
        {
            for (uint j=0; j<n_nodes; j++)
            {
                phi_xi_h(j,i) = phi_xih[j][i];
            }
        }
        
        // Shift quadrature points in -xi direction
        for (uint i=0; i<n_nodes; i++)
        {
            pts[i] = libMesh::Point(q_pts[i](0)-delta, q_pts[i](1), q_pts[i](2));
        }
        delete fe->_fe;
        fe->_fe = nullptr;
        fe->_initialized = false;
        fe->init(geom_elem, true, &pts);
        const std::vector<std::vector<Real>>& phi_xin = fe->get_phi();
        for (uint i=0; i<n_qps; i++)
        {
            for (uint j=0; j<n_nodes; j++)
            {
                phi_xi_n(j,i) = phi_xin[j][i];
            }
        }
        
        // Calculate second order central difference approximation to dphi_dxi
        RealMatrixX dphi_dxi_fd = RealMatrixX::Zero(n_nodes, n_qps);
        for (uint i=0; i<n_qps; i++) // Iterate Over Quadrature Points
        {
            for (uint j=0; j<n_nodes; j++) // Iterative Over Shape Functions
            {
                dphi_dxi_fd(j,i) = (phi_xi_h(j,i) - phi_xi_n(j,i))/(2.0*delta) ;
            }
        }
        //libMesh::out << "dphi_dxi:\n" << dphi_dxi_0 << std::endl;
        //libMesh::out << "dphi_dxi_fd:\n" << dphi_dxi_fd << std::endl;
        
        std::vector<double> dPhi_dxi =    eigen_matrix_to_std_vector(dphi_dxi_0);
        std::vector<double> dPhi_dxi_fd = eigen_matrix_to_std_vector(dphi_dxi_fd);
        
        REQUIRE_THAT( dPhi_dxi, Catch::Approx<double>(dPhi_dxi_fd) );
    }
    
    
    SECTION("shape function eta derivative finite difference check")
    {
        // Get the shape function derivative values at the quadrature points
        const std::vector<std::vector<libMesh::RealVectorValue>>& dphi = fe->get_dphi();
        
        uint n_qps = dphi[0].size();
        
        RealMatrixX dphi_deta_0 = RealMatrixX::Zero(n_nodes, n_qps);
        RealMatrixX phi_eta_h = RealMatrixX::Zero(n_nodes, n_qps);
        RealMatrixX phi_eta_n = RealMatrixX::Zero(n_nodes, n_qps);
        
        for (uint i=0; i<n_qps; i++)
        {
            for (uint j=0; j<n_nodes; j++)
            {
                dphi_deta_0(j,i) = dphi[j][i](1);
            }
        }
        
        // Get the quadrature points
        const std::vector<libMesh::Point>& q_pts = fe->get_qpoints();
        
        // Shift the Quadrature Points in the xi direction
        Real delta = 0.0001220703125; // sqrt(sqrt(eps))
        std::vector<libMesh::Point> pts;
        pts.reserve(n_nodes);
        
        // Shift quadrature points in +eta direction
        for (uint i=0; i<n_nodes; i++)
        {
            pts.push_back(libMesh::Point(q_pts[i](0), q_pts[i](1)+delta, q_pts[i](2)));
        }
        delete fe->_fe;
        fe->_fe = nullptr;
        fe->_initialized = false;
        fe->init(geom_elem, true, &pts);
        const std::vector<std::vector<Real>>& phi_etah = fe->get_phi();
        for (uint i=0; i<n_qps; i++)
        {
            for (uint j=0; j<n_nodes; j++)
            {
                phi_eta_h(j,i) = phi_etah[j][i];
            }
        }
        
        // Shift quadrature points in -eta direction
        for (uint i=0; i<n_nodes; i++)
        {
            pts[i] = libMesh::Point(q_pts[i](0), q_pts[i](1)-delta, q_pts[i](2));
        }
        delete fe->_fe;
        fe->_fe = nullptr;
        fe->_initialized = false;
        fe->init(geom_elem, true, &pts);
        const std::vector<std::vector<Real>>& phi_etan = fe->get_phi();
        for (uint i=0; i<n_qps; i++)
        {
            for (uint j=0; j<n_nodes; j++)
            {
                phi_eta_n(j,i) = phi_etan[j][i];
            }
        }
        
        // Calculate second order central difference approximation to dphi_dxi
        RealMatrixX dphi_deta_fd = RealMatrixX::Zero(n_nodes, n_qps);
        for (uint i=0; i<n_qps; i++) // Iterate Over Quadrature Points
        {
            for (uint j=0; j<n_nodes; j++) // Iterative Over Shape Functions
            {
                dphi_deta_fd(j,i) = (phi_eta_h(j,i) - phi_eta_n(j,i))/(2.0*delta) ;
            }
        }
        //libMesh::out << "dphi_deta:\n" << dphi_deta_0 << std::endl;
        //libMesh::out << "dphi_deta_fd:\n" << dphi_deta_fd << std::endl;
        
        std::vector<double> dPhi_deta =    eigen_matrix_to_std_vector(dphi_deta_0);
        std::vector<double> dPhi_deta_fd = eigen_matrix_to_std_vector(dphi_deta_fd);
        
        REQUIRE_THAT( dPhi_deta, Catch::Approx<double>(dPhi_deta_fd) );
    }
}
