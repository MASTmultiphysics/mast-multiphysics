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


// BOOST includes
#include <boost/test/unit_test.hpp>


// MAST includes
#include "examples/structural/bar_extension/bar_extension.h"
#include "tests/base/check_sensitivity.h"
#include "base/nonlinear_system.h"



// libMesh includes
#include "libmesh/numeric_vector.h"


BOOST_FIXTURE_TEST_SUITE  (Structural1DBarExtension,
                           MAST::BarExtension)

BOOST_AUTO_TEST_CASE   (BarExtensionSolution) {
    
    const Real
    tol      = 1.e-2;

    this->init(libMesh::EDGE2, false);
    this->solve();

    // check the solution
    // this is a constant stress extension problem with the
    // stress       =  pressure
    // strain       =  pressure/E
    // displacement =  strain * length
    //
    // iterate over each node, and compare the nodal solution with the
    // expected anlaytical value
    unsigned int
    dof_num = 0;
    libMesh::MeshBase::const_node_iterator
    it     =  _mesh->local_nodes_begin(),
    end    =  _mesh->local_nodes_end();
    
    Real
    press       = (*_press)(),
    Eval        = (*_E)(),
    analytical  = 0.,
    numerical   = 0.;
    
    for ( ; it!=end; it++) {
        const libMesh::Node* node = *it;
        dof_num      = node->dof_number(_sys->number(), _structural_sys->vars()[0],  0);
        analytical   =   (*node)(0) * press/Eval;
        numerical    =   _sys->solution->el(dof_num);
        BOOST_CHECK(MAST::compare_value(analytical, numerical, tol));
    }
    
    // make sure that each stress object has a single stored value
    for (unsigned int i=0; i<_outputs.size(); i++) {
        BOOST_CHECK(_outputs[i]->n_elem_in_storage() == 1);
    }
    
    // now check the stress value in each element, which should be the same as
    // the pressure value specified for the problem
    for (unsigned int i=0; i<_outputs.size(); i++) {
        // the call to all elements should actually include a single element only
        // the p-norm used is for p=2.
        numerical = _outputs[i]->von_Mises_p_norm_functional_for_all_elems(2);
        BOOST_CHECK(MAST::compare_value(press, numerical, tol));
    }
}


BOOST_AUTO_TEST_CASE   (BarExtensionSensitivity) {
    
    this->init(libMesh::EDGE2, false);
    MAST::check_sensitivity(*this);
}


BOOST_AUTO_TEST_SUITE_END()



