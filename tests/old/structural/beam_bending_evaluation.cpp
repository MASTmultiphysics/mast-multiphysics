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
#include "examples/structural/beam_bending/beam_bending.h"
#include "tests/base/check_sensitivity.h"
#include "base/nonlinear_system.h"


BOOST_FIXTURE_TEST_SUITE  (Structural1DBeamBending,
                           MAST::BeamBending)

BOOST_AUTO_TEST_CASE   (BeamBendingSolution) {
    
    const Real
    tol      = 1.e-2;

    this->init(libMesh::EDGE2, false);
    this->solve();
    
    // check the solution
    // iterate over each node, and compare the nodal solution with the
    // expected anlaytical value
    unsigned int
    dof_num = 0;
    libMesh::MeshBase::const_node_iterator
    it     =  _mesh->local_nodes_begin(),
    end    =  _mesh->local_nodes_end();
    
    Real
    press       = (*_press)(),
    th_y        = (*_thy)(),
    th_z        = (*_thz)(),
    Izz         = th_z*pow(th_y,3)/12.,
    x           = 0.,
    xi          = 0.,
    eta         = 0.,
    Eval        = (*_E)(),
    analytical  = 0.,
    numerical   = 0.;
    
    // analytical solution to the simply supported problem is
    // w(x)    = p/EI ( x^4/24 - L x^3/12 + L^2 x^2/24)
    // dwdx(x) = p/EI ( x^3/6  - L x^2/4  + L^2  x/12)
    
    for ( ; it!=end; it++) {
        const libMesh::Node* node = *it;
        x            = (*node)(0);

        // v-displacement
        analytical   = -press/Eval/Izz*(pow(x,4)/24. -
                                        _length*pow(x,3)/12. +
                                        pow(_length,2)*pow(x,2)/24.);
        
        dof_num      = node->dof_number(_sys->number(),
                                        _structural_sys->vars()[1], // v-displ.
                                        0);
        numerical    =   _sys->solution->el(dof_num);

        BOOST_CHECK(MAST::compare_value(analytical, numerical, tol));
        
        // theta-z rotation
        analytical   = -press/Eval/Izz*(pow(x,3)/6. -
                                        _length*pow(x,2)/4. +
                                        pow(_length,2)*x/12.);
        
        dof_num      = node->dof_number(_sys->number(),
                                        _structural_sys->vars()[5], // tz-rotation
                                        0);
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
        
        // get the element and the nodes to evaluate the stress
        const libMesh::Elem& e  = **(_outputs[i]->get_elem_subset().begin());
        
        const std::vector<MAST::StressStrainOutputBase::Data*>&
        data = _outputs[i]->get_stress_strain_data_for_elem(&e);
        
        // find the location of quadrature point
        for (unsigned int j=0; j<data.size(); j++) {

            // logitudinal strain for this location
            numerical = data[j]->stress()(0);
            
            xi   = data[j]->point_location_in_element_coordinate()(0);
            eta  = data[j]->point_location_in_element_coordinate()(1);
            
            // assuming linear Lagrange interpolation for elements
            x =  e.point(0)(0) * (1.-xi)/2. +  e.point(1)(0) * (1.+xi)/2.;
            // this gives the rotation at this point.
            analytical   = -press/Eval/Izz*(pow(x,2)/2. -
                                            _length*x/2. +
                                            pow(_length,2)/12.);
            // use this to evaluate the stress. The stress at both upper and
            // lower skins is same, but in opposite sign.
            analytical  *= -Eval*th_y*eta/2.;
            
            BOOST_CHECK(MAST::compare_value(analytical, numerical, tol));
        }
    }

}


BOOST_AUTO_TEST_CASE   (BeamBendingSensitivity) {
    
    
    this->init(libMesh::EDGE2, false);
    MAST::check_sensitivity(*this);
}


BOOST_AUTO_TEST_SUITE_END()

