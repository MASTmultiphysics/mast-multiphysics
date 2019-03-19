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
#include "examples/structural/membrane_extension_uniaxial_stress/membrane_extension_uniaxial.h"
#include "examples/structural/membrane_extension_biaxial_stress/membrane_extension_biaxial.h"
#include "tests/base/check_sensitivity.h"
#include "base/nonlinear_system.h"



void
uniaxial_membrane_sensitivity(MAST::MembraneExtensionUniaxial& mem) {

    const Real
    tol      = 1.e-2;

    // verify the sensitivity solution of this system
    RealVectorX
    sol,
    dsol,
    dsol_fd;
    
    mem.solve();
    
    // make sure that each stress object has a single stored value
    for (unsigned int i=0; i<mem._outputs.size(); i++)
        BOOST_CHECK((mem._outputs[i]->n_elem_in_storage() == 1));
    
    const Real
    press       = (*mem._press)(),
    Eval        = (*mem._E)(),
    nuval       = (*mem._nu)(),
    p_val       = 2.;
    
    Real
    analytical  = 0.,
    numerical   = 0.;

    unsigned int
    dof_num = 0;
    libMesh::MeshBase::const_node_iterator
    it     =  mem._mesh->local_nodes_begin(),
    end    =  mem._mesh->local_nodes_end();
    
    // now clear the stress data structures
    mem.clear_stresss();
    
    // now iterate over all the parameters and calculate the analytical
    // sensitivity and compare with the numerical sensitivity
    
    /////////////////////////////////////////////////////////
    // now evaluate the direct sensitivity
    /////////////////////////////////////////////////////////
    
    for (unsigned int i=0; i<mem._params_for_sensitivity.size(); i++ ) {
        
        MAST::Parameter& f = *mem._params_for_sensitivity[i];
        
        // calculate the analytical sensitivity
        // analysis is required at the baseline before sensitivity solution
        // and the solution has changed after the previous perturbed solution
        mem.solve();
        const libMesh::NumericVector<Real>& dsol_vec = mem.sensitivity_solve(f);
        
        // make sure that each stress object has a single stored value
        for (unsigned int i=0; i<mem._outputs.size(); i++)
            BOOST_CHECK((mem._outputs[i]->n_elem_in_storage() == 1));
        
        BOOST_TEST_MESSAGE("  ** dX/dp (total) wrt : " << f.name() << " **");
        it = mem._mesh->local_nodes_begin();
        for ( ; it!=end; it++) {
            const libMesh::Node* node = *it;
            
            // the x-displacement
            dof_num      = node->dof_number(mem._sys->number(),
                                            mem._structural_sys->vars()[0],
                                            0);
            analytical   = 0.;
            if (f.name() == "E")
                analytical   =   -(*node)(0) * press/pow(Eval,2);
            else if (f.name() == "nu")
                analytical   =   0.;
            else if (f.name() == "th")
                analytical   =   0.;
            else
                libmesh_error(); // should not get here
            numerical    =   dsol_vec.el(dof_num);
            BOOST_CHECK(MAST::compare_value(analytical, numerical, tol));
            
            // the y-displacement
            dof_num      = node->dof_number(mem._sys->number(),
                                            mem._structural_sys->vars()[1],
                                            0);
            analytical   = 0.;
            if (f.name() == "E")
                analytical   =   +nuval * (*node)(1) * press/pow(Eval,2);
            else if (f.name() == "nu")
                analytical   =   - (*node)(1) * press/Eval;
            else if (f.name() == "th")
                analytical   =   0.;
            else
                libmesh_error(); // should not get here
            numerical    =   dsol_vec.el(dof_num);
            BOOST_CHECK(MAST::compare_value(analytical, numerical, tol));
        }
        
        
        // now check the stress value in each element, which should be the same as
        // the pressure value specified for the problem
        BOOST_TEST_MESSAGE("  ** dvm-stress/dp (total) wrt : " << f.name() << " **");
        for (unsigned int i=0; i<mem._outputs.size(); i++) {
            // the call to all elements should actually include a single element only
            // the p-norm used is for p=2.
            numerical =
            mem._outputs[i]->von_Mises_p_norm_functional_sensitivity_for_all_elems(p_val, &f);
            BOOST_CHECK(MAST::compare_value(0., numerical, tol));
        }

        
        // now clear the stress data structures
        mem.clear_stresss();
    }
    
    
    /////////////////////////////////////////////////////////
    // now evaluate the adjoint sensitivity
    /////////////////////////////////////////////////////////
    
    

}



BOOST_FIXTURE_TEST_SUITE  (Structural2DMembraneExtensionUniaxial,
                           MAST::MembraneExtensionUniaxial)



BOOST_AUTO_TEST_CASE   (MembraneExtensionUniaxialSolution) {
    
    const Real
    tol      = 1.e-2;

    this->init(libMesh::QUAD4, false);
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
    nuval       = (*_nu)(),
    analytical  = 0.,
    numerical   = 0.;
    
    for ( ; it!=end; it++) {
        const libMesh::Node* node = *it;
        
        // the x-displacement
        dof_num      = node->dof_number(_sys->number(), _structural_sys->vars()[0],  0);
        analytical   =   (*node)(0) * press/Eval;
        numerical    =   _sys->solution->el(dof_num);
        BOOST_CHECK(MAST::compare_value(analytical, numerical, tol));

        // the y-displacement
        dof_num      = node->dof_number(_sys->number(), _structural_sys->vars()[1],  0);
        analytical   =   -nuval * (*node)(1) * press/Eval;
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


BOOST_AUTO_TEST_CASE   (MembraneExtensionUniaxialSensitivity) {
    
    uniaxial_membrane_sensitivity(*this);
}


BOOST_AUTO_TEST_SUITE_END()



BOOST_FIXTURE_TEST_SUITE  (Structural2DMembraneExtensionBiaxial,
                           MAST::MembraneExtensionBiaxial)

BOOST_AUTO_TEST_CASE   (MembraneExtensionBiaxialSensitivity) {
    
    MAST::check_sensitivity(*this);
}


BOOST_AUTO_TEST_SUITE_END()

