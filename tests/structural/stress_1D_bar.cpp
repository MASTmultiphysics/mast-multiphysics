/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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
#include "tests/structural/build_structural_elem_1D.h"
#include "tests/base/test_comparisons.h"
#include "tests/structural/build_structural_elem_1D.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/stress_output_base.h"


// libMesh includes
#include "libmesh/dof_map.h"


BOOST_FIXTURE_TEST_SUITE  (StructuralStressEvaluation, MAST::BuildStructural1DElem)

BOOST_AUTO_TEST_CASE   (Stress1DBar) {
    

    
    // stress output
    MAST::StressStrainOutputBase output;
    std::multimap<libMesh::subdomain_id_type, MAST::OutputFunctionBase*> output_map;
    output_map.insert(std::pair<libMesh::subdomain_id_type, MAST::OutputFunctionBase*>
                      (0, &output));

    
    // get reference to the element in this mesh
    const libMesh::Elem& elem = **(_mesh->local_elements_begin());
    
    // now create the structural element
    std::auto_ptr<MAST::StructuralElementBase>
    e(MAST::build_structural_element(*_structural_sys, elem, *_p_card, false).release());
    
    
    // number of dofs in this element
    const libMesh::DofMap& dofmap = _sys->get_dof_map();
    std::vector<unsigned int> dof_ids;
    dofmap.dof_indices(&elem, dof_ids);
    
    const unsigned int ndofs = (unsigned int)dof_ids.size();
    
    // now get the residual and Jacobian evaluations
    RealVectorX
    x0           = RealVectorX::Zero(ndofs),
    xdot0        = RealVectorX::Zero(ndofs),
    x            = RealVectorX::Zero(ndofs),
    xdot         = RealVectorX::Zero(ndofs),
    stress0      = RealVectorX::Zero(6),
    stress       = RealVectorX::Zero(6),
    strain0      = RealVectorX::Zero(6),
    strain       = RealVectorX::Zero(6),
    dvm_dX0      = RealVectorX::Zero(ndofs),
    dvm_dX_fd    = RealVectorX::Zero(ndofs),
    dstressdp    = RealVectorX::Zero(6),
    dstraindp    = RealVectorX::Zero(6);
    
    RealMatrixX
    dstressdX0   = RealMatrixX::Zero(6, ndofs),
    dstressdX_fd = RealMatrixX::Zero(6, ndofs),
    dstraindX0   = RealMatrixX::Zero(6, ndofs),
    dstraindX_fd = RealMatrixX::Zero(6, ndofs);
    
    const Real
    delta   = 1.0e-5,
    tol     = 1.0e-7;
    
    Real
    vm0     = 0.,
    vm      = 0.,
    dvmdp   = 0.,
    p0      = 0.,
    dp      = 0.;
    
    // tell the element about the solution and velocity
    x0(1) = 2.e-2;          // some arbitrary axial deformation
    e->set_solution(x0);
    e->set_velocity(xdot0);
    
    // also ask for the derivative of the stress/strain
    e->volume_output_quantity(true, false, output_map);

    // get access to the vector of stress/strain data for this element.
    {
        const std::vector<MAST::StressStrainOutputBase::Data*>&
        stress_data = output.get_stress_strain_data();
        
        libmesh_assert_equal_to(stress_data.size(), 1); // this should have one element
        
        stress0     = stress_data[0]->stress();
        strain0     = stress_data[0]->strain();
        dstressdX0  = stress_data[0]->get_dstress_dX();
        dstraindX0  = stress_data[0]->get_dstrain_dX();
        vm0         = stress_data[0]->von_Mises_stress();
        dvm_dX0     = stress_data[0]->dvon_Mises_stress_dX();
        
        output.clear();
    }
    
    // now, verify the adjoints using finite differencing
    for (unsigned int i=0; i<ndofs; i++) {
        x      = x0;
        x(i)  += delta;
        e->set_solution(x);
        e->set_velocity(xdot);
        
        e->volume_output_quantity(false, false, output_map);

        // now use the updated stress to calculate the finite difference data
        {
            const std::vector<MAST::StressStrainOutputBase::Data*>&
            stress_data = output.get_stress_strain_data();
            
            libmesh_assert_equal_to(stress_data.size(), 1); // this should have one element
            
            stress              = stress_data[0]->stress();
            strain              = stress_data[0]->strain();
            dstressdX_fd.col(i) = (stress-stress0)/delta;
            dstraindX_fd.col(i) = (strain-strain0)/delta;
            vm                  = stress_data[0]->von_Mises_stress();
            dvm_dX_fd(i)        = (vm-vm0)/delta;
            
            output.clear();
        }
    }
    
    // now compare the matrices
    BOOST_CHECK(MAST::compare_matrix(   dstressdX0,    dstressdX_fd, tol));
    BOOST_CHECK(MAST::compare_matrix(   dstraindX0,    dstraindX_fd, tol));
    BOOST_CHECK(MAST::compare_vector(      dvm_dX0,       dvm_dX_fd, tol));
    
    
    // now, check the sensitivity with respect to various parameters
    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    std::vector<MAST::Parameter*> params_for_sensitivity;
    params_for_sensitivity.push_back(_E);
    params_for_sensitivity.push_back(_nu);
    params_for_sensitivity.push_back(_thy);
    params_for_sensitivity.push_back(_thz);
    
    for (unsigned int i=0; i<params_for_sensitivity.size(); i++) {
        
        MAST::Parameter& f = *params_for_sensitivity[i];
        
        // set the sensitivity of solution to be zero
        x     = RealVectorX::Zero(ndofs);
        e->sensitivity_param  = &f;
        e->set_solution(x0);
        e->set_solution(x, true);  // sets the sensitivity to zero
        
        // first check the partial derivative of the quantity wrt the parameter
        p0           = f();
        if (p0 > 0)
            dp = delta * p0;
        else
            dp = delta;
        f()         += dp;
        
        // get the stress, strain and von Mises stress sensitivity
        e->volume_output_quantity(false, true, output_map);

        // reset the parameter value
        f()        = p0;
        
        // next, check the total derivative of the quantity wrt the parameter
        {
            const std::vector<MAST::StressStrainOutputBase::Data*>&
            stress_data = output.get_stress_strain_data();
            
            libmesh_assert_equal_to(stress_data.size(), 1); // this should have one element
            
            dstressdp           = stress_data[0]->get_stress_sensitivity(&f);
            dstraindp           = stress_data[0]->get_strain_sensitivity(&f);
            dvmdp               = stress_data[0]->dvon_Mises_stress_dp  (&f);
            

            stress              = (stress_data[0]->stress() - stress0)/dp;
            strain              = (stress_data[0]->strain() - strain0)/dp;
            vm                  = (stress_data[0]->von_Mises_stress() - vm0)/dp;

            output.clear();
        }
        
        // compare and benchmark the values
        BOOST_CHECK(MAST::compare_vector(dstressdp, stress, tol));
        BOOST_CHECK(MAST::compare_vector(dstraindp, strain, tol));
        BOOST_CHECK(MAST::compare_value (    dvmdp,     vm, tol));
    }
    
}

BOOST_AUTO_TEST_SUITE_END()


