#include "catch.hpp"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/distributed_mesh.h"
#include "libmesh/face_quad4.h"
#include "libmesh/edge_edge2.h"

// Custom includes
#include "test_helpers.h"

extern libMesh::LibMeshInit* p_global_init;

TEST_CASE("libmesh_mesh_generation_1d",
          "[mesh][1D]")
{
    const int n_elems = 1;
    const int n_nodes = 2;
    const int n_dim = 1;
    
    RealMatrixX temp = RealMatrixX::Zero(3,n_nodes);
    temp << -1.0, 1.0, 0.0, 0.0, 0.0, 0.0;
    const RealMatrixX X = temp;
    
    SECTION("replicated_mesh_1d")
    {
        /**
        *  First create the mesh with the one element we are testing.
        */
        // Setup the mesh properties
        libMesh::ReplicatedMesh mesh(p_global_init->comm());
        mesh.set_mesh_dimension(n_dim);
        mesh.set_spatial_dimension(n_dim);
        mesh.reserve_elem(n_elems);
        mesh.reserve_nodes(n_nodes);
        
        // Add nodes to the mesh
        for (uint i=0; i<n_nodes; i++)
        {
            mesh.add_point(libMesh::Point(X(0,i), X(1,i), X(2,i)));
        }
        
        // Add the element to the mesh
        libMesh::Elem *reference_elem = new libMesh::Edge2;
        reference_elem->set_id(0);    
        reference_elem->subdomain_id() = 0;
        reference_elem = mesh.add_elem(reference_elem);
        for (int i=0; i<n_nodes; i++)
        {
            reference_elem->set_node(i) = mesh.node_ptr(i);
        }
        
        // Prepare the mesh for use
        mesh.prepare_for_use();
        
        // Get the volume calculated by libMesh (actually length for 1D elements)
        const Real elem_volume = reference_elem->volume();
        
        // Set the true length
        const Real true_volume = 2.0;
        
        // Ensure the libMesh element has the expected volume
        REQUIRE( elem_volume == true_volume );
    }
}


TEST_CASE("libmesh_mesh_generation_2d",
          "[mesh],[2D]")
{
    const int n_elems = 1;
    const int n_nodes = 4;
    const int n_dim = 2;
    
    // Point Coordinates
    RealMatrixX temp = RealMatrixX::Zero(3,4);
    temp << -1.0,  1.0, 1.0, -1.0, 
            -1.0, -1.0, 1.0,  1.0, 
             0.0,  0.0, 0.0,  0.0;
    const RealMatrixX X = temp;
    
    SECTION("replicated_mesh")
    {
        /**
        *  First create the mesh with the one element we are testing.
        */
        // Setup the mesh properties
        libMesh::ReplicatedMesh mesh(p_global_init->comm());
        mesh.set_mesh_dimension(n_dim);
        mesh.set_spatial_dimension(n_dim);
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
        
        // Get the volume calculated by libMesh (actually area for 2D elements)
        const Real elem_volume = reference_elem->volume();
        
        // Calculate true area using 2D shoelace formula
        const Real true_volume = get_shoelace_area(X);
        
        // Ensure the libMesh element has the expected volume
        REQUIRE( elem_volume == true_volume );
    }
    
   // TODO Need to test distributed mesh generation as well
}
