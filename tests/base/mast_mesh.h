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

#ifndef __test__mast_mesh__
#define __test__mast_mesh__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/distributed_mesh.h"
#include "libmesh/face_quad4.h"
#include "libmesh/edge_edge2.h"

// Test includes
#include "catch.hpp"
#include "test_helpers.h"

extern libMesh::LibMeshInit* p_global_init;

namespace TEST {

    /**
     * Storage class for a mesh consisting of a single element used in testing.
     *
     * The single element has an ID of 0 (zero) and is placed into the subdomain of ID 0 (zero).
     */
    class TestMeshSingleElement {
    public:
        int n_elems; ///< Number of elements in the test mesh
        int n_nodes; ///< Number of nodes per element in the test mesh
        int n_dim;   ///< Dimension of the test element (1, 2, 3)
        libMesh::Elem* reference_elem; ///< Pointer to the actual libMesh element object
        libMesh::ReplicatedMesh mesh;  ///< The actual libMesh mesh object
        // Convert this to pointer to enable both Replicated/Distributed Mesh
        // libMesh::UnstructuredMesh* mesh;
        // ---> currently can't run this with DistributedMesh. On second processor this.reference_elem doesn't exist!

        /**
         * Construct a single element mesh using the specified type and nodal coordinates. Nodal
         * connectivity is specified according to the Exodus-II mesh format. Valid element types
         * are currently:
         *  - 1D: EDGE2
         *  - 2D: QUAD4
         *
         * @param e_type libMesh element type to create
         * @param coordinates (3 by n) matrix where each column specifies a node with rows giving
         *                    the x, y, z locations
         */
        TestMeshSingleElement(libMesh::ElemType e_type, RealMatrixX& coordinates):
            mesh(p_global_init->comm()) {
            n_elems = 1;
            n_nodes = coordinates.cols();

            mesh.reserve_elem(n_elems);
            mesh.reserve_nodes(n_nodes);
            mesh.set_spatial_dimension(3);

            for (auto i = 0; i < n_nodes; i++) {
                mesh.add_point(libMesh::Point(coordinates(0,i),coordinates(1,i), coordinates(2,i)), i, 0);
            }

            switch (e_type) {
                case libMesh::EDGE2:
                    n_dim = 1;
                    reference_elem = new libMesh::Edge2;
                    break;
                case libMesh::QUAD4:
                    n_dim = 2;
                    reference_elem = new libMesh::Quad4;
                    break;
                default:
                    libmesh_error_msg("Invalid element type; " << __PRETTY_FUNCTION__
                                                               << " in " << __FILE__ << " at line number " << __LINE__);
            }

            mesh.set_mesh_dimension(n_dim);
            reference_elem->set_id(0);
            reference_elem->subdomain_id() = 0;
            reference_elem = mesh.add_elem(reference_elem);

            for (int i=0; i<n_nodes; i++) {
                reference_elem->set_node(i) = mesh.node_ptr(i);
            }

            mesh.prepare_for_use();
        };

        /**
         * Update the nodal coordinates in the mesh.
         *
         * @param new_coordinates (3 by n) matrix where each column specifies a node with rows
         *                        giving the x, y, z locations
         */
        void update_coordinates(RealMatrixX& new_coordinates) {
            for (int i=0; i<n_nodes; i++)
            {
                *mesh.node_ptr(i) = libMesh::Point(new_coordinates(0,i), new_coordinates(1,i), new_coordinates(2,i));
            }
        };
    };
} // TEST namespace

#endif // __test__mast_mesh__