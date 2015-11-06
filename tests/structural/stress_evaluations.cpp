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
#include "tests/structural/build_structural_elem_2D.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "tests/base/test_comparisons.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/stress_output_base.h"


// libMesh includes
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"

extern const Real
delta,
tol;

void
set_deformation(const unsigned int dim,
                RealVectorX& vec) {
    if (dim == 1)
        // assuming an edge2 with Lagrange interpolation
        vec(1) = 2.e-2;          // some arbitrary axial deformation
    else if (dim == 2) {
        
        // assuming a quad4 with Lagrange interpolation
        vec(2) = 2.e-2;
        vec(3) = 4.e-2;
        
        vec(6) = -2.e-2;
        vec(7) = +6.e-2;
    }

}

template <typename ValType>
void check_stress (ValType& v) {

    // stress output
    MAST::StressStrainOutputBase output;
    std::multimap<libMesh::subdomain_id_type, MAST::OutputFunctionBase*> output_map;
    output_map.insert(std::pair<libMesh::subdomain_id_type, MAST::OutputFunctionBase*>
                      (0, &output));
    
    
    // get reference to the element in this mesh
    const libMesh::Elem& elem = **(v._mesh->local_elements_begin());
    
    // now create the structural element
    std::auto_ptr<MAST::StructuralElementBase>
    e(MAST::build_structural_element(*v._structural_sys,
                                     elem,
                                     *v._p_card).release());
    
    
    // number of dofs in this element
    const libMesh::DofMap& dofmap = v._sys->get_dof_map();
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
    dvmf_dX0     = RealVectorX::Zero(ndofs),
    dvm_dX_fd    = RealVectorX::Zero(ndofs),
    dvmf_dX_fd   = RealVectorX::Zero(ndofs),
    dstressdp    = RealVectorX::Zero(6),
    dstraindp    = RealVectorX::Zero(6);
    
    RealMatrixX
    dstressdX0   = RealMatrixX::Zero(6, ndofs),
    dstressdX_fd = RealMatrixX::Zero(6, ndofs),
    dstraindX0   = RealMatrixX::Zero(6, ndofs),
    dstraindX_fd = RealMatrixX::Zero(6, ndofs);
    
    
    Real
    vm0          = 0.,
    vm           = 0.,
    vmf0         = 0.,
    dvmdp        = 0.,
    dvmf_dp      = 0.,
    dvmf_dp_fd   = 0.,
    p0           = 0.,
    dp           = 0.;
    
    const Real
    pval    = 3.;
    
    // tell the element about the solution and velocity
    set_deformation(elem.dim(), x0);
    e->set_solution(x0);
    e->set_velocity(xdot0);
    
    // also ask for the derivative of the stress/strain
    e->volume_output_quantity(true, false, output_map);
    
    // get access to the vector of stress/strain data for this element.
    {
        const std::vector<MAST::StressStrainOutputBase::Data*>&
        stress_data = output.get_stress_strain_data_for_elem(&elem);
        
        libmesh_assert_equal_to(stress_data.size(), 1); // this should have one element
        
        stress0     = stress_data[0]->stress();
        strain0     = stress_data[0]->strain();
        dstressdX0  = stress_data[0]->get_dstress_dX();
        dstraindX0  = stress_data[0]->get_dstrain_dX();
        vm0         = stress_data[0]->von_Mises_stress();
        dvm_dX0     = stress_data[0]->dvon_Mises_stress_dX();
        vmf0        = output.von_Mises_p_norm_functional_for_all_elems(pval);
        dvmf_dX0    = output.von_Mises_p_norm_functional_state_derivartive_for_all_elems(pval);
        
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
            stress_data = output.get_stress_strain_data_for_elem(&elem);
            
            libmesh_assert_equal_to(stress_data.size(), 1); // this should have one element
            
            stress              = stress_data[0]->stress();
            strain              = stress_data[0]->strain();
            dstressdX_fd.col(i) = (stress-stress0)/delta;
            dstraindX_fd.col(i) = (strain-strain0)/delta;
            vm                  = stress_data[0]->von_Mises_stress();
            dvm_dX_fd(i)        = (vm-vm0)/delta;
            dvmf_dX_fd(i)       = (output.von_Mises_p_norm_functional_for_all_elems(pval)-vm0)/delta;
            
            output.clear();
        }
    }
    
    // now compare the matrices
    BOOST_TEST_MESSAGE("  ** Stress Derivative  wrt State **");
    BOOST_CHECK(MAST::compare_matrix(   dstressdX0,    dstressdX_fd, tol));
    BOOST_TEST_MESSAGE("  ** Strain Derivative  wrt State **");
    BOOST_CHECK(MAST::compare_matrix(   dstraindX0,    dstraindX_fd, tol));
    BOOST_TEST_MESSAGE("  ** VM-Stress Derivative  wrt State **");
    BOOST_CHECK(MAST::compare_vector(      dvm_dX0,       dvm_dX_fd, tol));
    BOOST_TEST_MESSAGE("  ** VM_func-Stress Derivative  wrt State **");
    BOOST_CHECK(MAST::compare_vector(      dvmf_dX0,       dvmf_dX_fd, tol));
    
    // now, check the sensitivity with respect to various parameters
    
    // NOTE:
    // First, we will set the dX/dp = 0, so that the derivative of the output
    // quantity is only the partial derivative wrt the parameter.
    //
    
    for (unsigned int i=0; i<v._params_for_sensitivity.size(); i++) {
        
        MAST::Parameter& f = *v._params_for_sensitivity[i];
        
        // set the sensitivity of solution to be zero
        x     = RealVectorX::Zero(ndofs);
        e->sensitivity_param  = &f;
        e->set_solution(x0);
        e->set_solution(x, true);  // sets the sensitivity to zero
        
        // identify the perturbation in the parameter
        p0           = f();
        (p0 > 0)?  dp=delta*p0 : dp=delta;
        f()         += dp;
        
        // get the stress, strain and von Mises stress sensitivity
        e->volume_output_quantity(false, true, output_map);
        
        // reset the parameter value
        f()        = p0;
        
        // next, check the total derivative of the quantity wrt the parameter
        {
            const std::vector<MAST::StressStrainOutputBase::Data*>&
            stress_data = output.get_stress_strain_data_for_elem(&elem);
            
            libmesh_assert_equal_to(stress_data.size(), 1); // this should have one element
            
            dstressdp           = stress_data[0]->get_stress_sensitivity(&f);
            dstraindp           = stress_data[0]->get_strain_sensitivity(&f);
            dvmdp               = stress_data[0]->dvon_Mises_stress_dp  (&f);
            dvmf_dp             =
            output.von_Mises_p_norm_functional_sensitivity_for_all_elems(pval, &f);
            
            stress              = (stress_data[0]->stress() - stress0)/dp;
            strain              = (stress_data[0]->strain() - strain0)/dp;
            vm                  = (stress_data[0]->von_Mises_stress() - vm0)/dp;
            dvmf_dp_fd          =
            (output.von_Mises_p_norm_functional_for_all_elems(pval)-vmf0)/dp;
            
            output.clear();
        }
        
        // compare and benchmark the values
        BOOST_TEST_MESSAGE("  ** dStress/dp (partial) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_vector(dstressdp, stress, tol));
        BOOST_TEST_MESSAGE("  ** dStrain/dp (partial) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_vector(dstraindp, strain, tol));
        BOOST_TEST_MESSAGE("  ** dVM-Stress/dp (partial) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_value (    dvmdp,     vm, tol));
        BOOST_TEST_MESSAGE("  ** dVM_func-Stress/dp (partial) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_value (    dvmf_dp,     dvmf_dp_fd, tol));
    }
    
    
    
    
    // NOTE:
    // Next, we will set the dp = 0, so that the derivative of the output
    // quantity is only the partial derivative including sensitivity of
    // the state vector wrt the parameter.
    //
    
    for (unsigned int i=0; i<v._params_for_sensitivity.size(); i++) {
        
        MAST::Parameter& f = *v._params_for_sensitivity[i];
        
        // identify the perturbation in the parameter
        // if delta=1.e-5, and dp=1.e-3, then delta/dp = 1.e-2
        // magnitude of delta/dp should not matter, since the
        // sensitivity is linear in it. We really want delta and dp to be
        // small so that finite differencing can be accurately done
        p0           = f();
        (p0 > 0)?  dp=std::max(delta*p0, delta) : dp=delta;
        f()         += dp;
        
        // set the sensitivity of solution to be zero
        e->sensitivity_param  = &f;
        
        // now perturb each element of the vector independently
        for (unsigned int j=0; j<ndofs; j++) {
            
            // first perturb the solution for fd sensitivity
            x       = x0;
            x(j)   += delta;
            e->set_solution(x);
            // now, set the solution sensitivity for analytical sensitivity
            x       = RealVectorX::Zero(ndofs);
            x(j)    = delta/dp;
            e->set_solution(x, true);  // sets the sensitivity to zero
            
            
            
            // get the stress, strain and von Mises stress sensitivity
            e->volume_output_quantity(false, true, output_map);
            
            // next, check the total derivative of the quantity wrt the parameter
            {
                const std::vector<MAST::StressStrainOutputBase::Data*>&
                stress_data = output.get_stress_strain_data_for_elem(&elem);
                
                libmesh_assert_equal_to(stress_data.size(), 1); // this should have one element
                
                dstressdp           = stress_data[0]->get_stress_sensitivity(&f);
                dstraindp           = stress_data[0]->get_strain_sensitivity(&f);
                dvmdp               = stress_data[0]->dvon_Mises_stress_dp  (&f);
                dvmf_dp             =
                output.von_Mises_p_norm_functional_sensitivity_for_all_elems(pval, &f);

                
                stress              = (stress_data[0]->stress() - stress0)/dp;
                strain              = (stress_data[0]->strain() - strain0)/dp;
                vm                  = (stress_data[0]->von_Mises_stress() - vm0)/dp;
                dvmf_dp_fd          =
                (output.von_Mises_p_norm_functional_for_all_elems(pval)-vmf0)/dp;
                
                output.clear();
            }
            
            // compare and benchmark the values
            BOOST_TEST_MESSAGE("  ** dStress/dp (total) wrt : " << f.name() << " **");
            BOOST_CHECK(MAST::compare_vector(dstressdp, stress, tol));
            BOOST_TEST_MESSAGE("  ** dStrain/dp (total) wrt : " << f.name() << " **");
            BOOST_CHECK(MAST::compare_vector(dstraindp, strain, tol));
            BOOST_TEST_MESSAGE("  ** dVM-Stress/dp (total) wrt : " << f.name() << " **");
            BOOST_CHECK(MAST::compare_value (    dvmdp,     vm, tol));
            BOOST_TEST_MESSAGE("  ** dVM_func-Stress/dp (partial) wrt : " << f.name() << " **");
            BOOST_CHECK(MAST::compare_value (    dvmf_dp,     dvmf_dp_fd, tol));
        }
        
        // reset the parameter value
        f()        = p0;
    }
}


BOOST_FIXTURE_TEST_SUITE  (Structural1DStressEvaluation, MAST::BuildStructural1DElem)

BOOST_AUTO_TEST_CASE   (Stress1D) {
    
    check_stress(*this);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE  (Structural2DStressEvaluation, MAST::BuildStructural2DElem)

BOOST_AUTO_TEST_CASE   (Stress2D) {
    
    check_stress(*this);
}

BOOST_AUTO_TEST_SUITE_END()


