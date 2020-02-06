/** 
 * define below is needed to be able to access the 
 * initialize_green_lagrange_strain_operator which is protected. Be careful
 * though as the line belong can cause unexpected issues.
 */
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
#include "numerics/fem_operator_matrix.h"



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
TEST_CASE("quad4_linear_structural_strain_displacement_matrix", 
          "[quad],[quad4],[linear],][structural],[2D],[element]")
{
    const int n_elems = 1;
    const int n_nodes = 4;
    const uint n_dofs_per_node = 6;
    
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
    
    MAST::FEMOperatorMatrix
        Bmat_lin,
        Bmat_nl_x,
        Bmat_nl_y,
        Bmat_nl_u,
        Bmat_nl_v,
        Bmat_bend,
        Bmat_vk;
    
    const uint n_phi = (unsigned int)fe->get_phi().size();
    const uint n1 = elem->n_direct_strain_components();
    const uint n2 = 6*n_phi;
    const uint n3 = elem->n_von_karman_strain_components();
    
    RealMatrixX mat_x = RealMatrixX::Zero(3,2);
    RealMatrixX mat_y = RealMatrixX::Zero(3,2);
    
    RealVectorX strain = RealVectorX::Zero(3);
    
    uint qp = 0;
    
    /**
     * elem->initialize_green_lagrange_strain_operator method populates the
     * Bmat_lin, Bmat_nl_x, Bmat_nl_y, Bmat_nl_u, and Bmat_nl_v matrices.
     */
    Bmat_lin.reinit(n1, structural_system.n_vars(), n_phi); // three stress-strain components
    Bmat_nl_x.reinit(2, structural_system.n_vars(), n_phi);
    Bmat_nl_y.reinit(2, structural_system.n_vars(), n_phi);
    Bmat_nl_u.reinit(2, structural_system.n_vars(), n_phi);
    Bmat_nl_v.reinit(2, structural_system.n_vars(), n_phi);
    
    elem->initialize_green_lagrange_strain_operator(qp, *fe, elem_solution,
        strain, mat_x, mat_y, Bmat_lin, Bmat_nl_x, Bmat_nl_y, Bmat_nl_u,
        Bmat_nl_v);
    
    /**
     * std::unique_ptr<MAST::BendingOperator2D> bend;
     * bend->initialize_bending_strain_operator method populates the Bmat_bend
     * matrix. This part only exists if bending exists in the model.
     */
    Bmat_bend.reinit(n1, structural_system.n_vars(), n_phi);
    
    /**
     * elem->initialize_von_karman_strain_operator method populates the Bmat_vk
     * matrix. This part only exists if bending exists in the model AND 
     * nonlinear strains exist in the model.
     */
    Bmat_vk.reinit(n3, structural_system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    
    const std::vector<std::vector<libMesh::RealVectorValue>>& dphi = fe->get_dphi();
    
    SECTION("Linear in-plane strain displacement matrix size")
    {
        REQUIRE( Bmat_lin.m() == 3 ); // three strains, e_xx, e_yy, e_xy
        REQUIRE( Bmat_lin.n() == n_nodes * n_dofs_per_node);
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
        for (uint i=0; i<n_nodes; i++)
        {
            Bmat_lin_true(0,i)          = dphi[i][qp](0);
            Bmat_lin_true(2,i+n_nodes)  = dphi[i][qp](0);
            
            Bmat_lin_true(1,i+n_nodes)  = dphi[i][qp](1);
            Bmat_lin_true(2,i)          = dphi[i][qp](1);
        }
                
        std::vector<double> Bmat_lin_test =    eigen_matrix_to_std_vector(Bmat_lin_mat);
        std::vector<double> Bmat_lin_required = eigen_matrix_to_std_vector(Bmat_lin_true);
        
        REQUIRE_THAT( Bmat_lin_test, Catch::Approx<double>(Bmat_lin_required) );
    }
}
