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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

// MAST includes
#include "base/mast_data_types.h"
#include "catch.hpp"

// Test includes
#include "test_helpers.h"
#include "base/mast_mesh.h"


TEST_CASE("libmesh_mesh_generation_1d",
          "[mesh][1D]")
{
    RealMatrixX coords = RealMatrixX::Zero(3, 2);
    coords << -1.0, 1.0, 0.0,
               0.0, 0.0, 0.0;

    SECTION("Creation of a single 2-node 1D element (EDGE2) in a libMesh::ReplicatedMesh")
    {
        TEST::TestMeshSingleElement test_mesh(libMesh::EDGE2, coords);
        // Compare element against true volume (actually length for 1D elements)
        REQUIRE(test_mesh.reference_elem->volume() == 2.0);
    }
}


TEST_CASE("libmesh_mesh_generation_2d",
          "[mesh],[2D]")
{
    RealMatrixX coords = RealMatrixX::Zero(3,4);
    coords << -1.0,  1.0, 1.0, -1.0,
              -1.0, -1.0, 1.0,  1.0,
               0.0,  0.0, 0.0,  0.0;

    SECTION("Creation of a single 4-node 2D element (QUAD4) in a libMesh::ReplicatedMesh")
    {
        TEST::TestMeshSingleElement test_mesh(libMesh::QUAD4, coords);
        // Compare element against true volume (actually area for 2D elements)
        const Real true_volume = TEST::get_shoelace_area(coords);
        REQUIRE(test_mesh.reference_elem->volume() == true_volume);
    }
    
   // TODO Need to test distributed mesh generation as well.
   //   Note this will require some special consideration to account for actual
   //   parallel nature of the mesh.
}
