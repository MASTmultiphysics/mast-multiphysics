/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MAST_TESTS
#include <boost/test/unit_test.hpp>


// C++ includes
#include <memory>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/libmesh.h"


libMesh::LibMeshInit     *_mast_init         = nullptr;

struct GlobalTestFixture {
    
    GlobalTestFixture() {
        
        // create the libMeshInit function
        _mast_init = new libMesh::LibMeshInit(boost::unit_test::framework::master_test_suite().argc,
                                              boost::unit_test::framework::master_test_suite().argv);
    }
    
    ~GlobalTestFixture() {
        
        delete _mast_init;
    }
    
};


BOOST_GLOBAL_FIXTURE( GlobalTestFixture );



