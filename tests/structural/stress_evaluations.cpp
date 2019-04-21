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
#include "base/nonlinear_system.h"


// libMesh includes
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/edge_edge2.h"



void
set_deformation(const unsigned int dim,
                const unsigned int case_num,
                const libMesh::ElemType e,
                RealVectorX& vec) {
    if (dim == 1) {
        
        // assuming an edge2 with Lagrange interpolation
        vec   = RealVectorX::Zero(12);
        
        
        switch (case_num) {
            case 0: // some arbitrary axial deformation
                vec(1)  = 2.e-2;
                break;
            case 1: // some arbitrary bending deformation
                vec(9)  = 2.e-2;
                break;
            case 2: { // combination of axial and bending
                vec(1)  = 2.e-2;
                vec(9)  = 2.e-2;
            }
            case 3: { // combination of axial and bending
                vec(1)  = 2.e-1;   // u
                vec(3)  = 2.e-1;   // v
                vec(5)  = 2.e-1;   // w
                vec(9)  = 2.e-1;   // theta-y
                vec(11) = 2.e-1;   // theta-z
            }
                break;
            default:
                libmesh_error(); // should not get here
        }
    }
    else if (dim == 2) {
        
        // assuming a quad4 with Lagrange interpolation
        
        std::vector<unsigned int>
        mem_dofs,
        bend_dofs;
        
        if (e == libMesh::QUAD4) {
            vec   = RealVectorX::Zero(24);
            mem_dofs  = {2, 3, 6, 7};
            bend_dofs = {10, 11, 14, 15, 18, 19};
        }
        else if (e == libMesh::TRI3) {
            vec   = RealVectorX::Zero(18);
            mem_dofs  = {1, 2, 4, 5};
            bend_dofs = {7, 8, 10, 11, 13, 14};
        }
        else
            libmesh_error(); // should not get here.
        
        switch (case_num) {
                
            case 0: { // sets values for u and v
                vec(mem_dofs[0]) = 2.e-2;
                vec(mem_dofs[1]) = 4.e-2;
                
                vec(mem_dofs[2]) = -2.e-2;
                vec(mem_dofs[3]) = +6.e-2;
            }
                break;

            case 1: { // sets values for theta-x and theta-y
                vec(bend_dofs[0]) = 2.e-2;
                vec(bend_dofs[1]) = 4.e-2;
                
                vec(bend_dofs[2]) = -2.e-2;
                vec(bend_dofs[3]) = +6.e-2;

                vec(bend_dofs[4]) = -2.e-2;
                vec(bend_dofs[5]) = +6.e-2;
            }
                break;

            case 2: {
                vec(mem_dofs[0]) = 2.e-2;
                vec(mem_dofs[1]) = 4.e-2;
                
                vec(mem_dofs[2]) = -2.e-2;
                vec(mem_dofs[3]) = +6.e-2;

                vec(bend_dofs[0]) = 2.e-2;
                vec(bend_dofs[1]) = 4.e-2;
                
                vec(bend_dofs[2]) = -2.e-2;
                vec(bend_dofs[3]) = +6.e-2;

                vec(bend_dofs[4]) = -2.e-2;
                vec(bend_dofs[5]) = +6.e-2;
            }
                break;

            case 3: {
                vec(mem_dofs[0]) = 2.e-1;
                vec(mem_dofs[1]) = 4.e-1;
                
                vec(mem_dofs[2]) = -2.e-1;
                vec(mem_dofs[3]) = +6.e-1;
                
                vec(bend_dofs[0]) = 2.e-1;
                vec(bend_dofs[1]) = 4.e-1;
                
                vec(bend_dofs[2]) = -2.e-1;
                vec(bend_dofs[3]) = +6.e-1;
                
                vec(bend_dofs[4]) = -2.e-1;
                vec(bend_dofs[5]) = +6.e-1;
            }
                break;
                
            default:
                libmesh_assert(false); // should not get here.
        }
    }

}

template <typename ValType>
void check_stress (ValType& v, const RealVectorX& x0) {

    const Real
    delta    = 1.e-5,
    tol      = 1.e-2;

    // stress output
    MAST::StressStrainOutputBase output;
    std::vector<libMesh::Point> pts(1);
    pts[0] = libMesh::Point(0., 1., 0.);  // upper surface, elem mid-point

    //pts[1] = libMesh::Point(0.,-1., 0.);  // lower surface, elem mid-point
    output.set_points_for_evaluation(pts);
    std::multimap<libMesh::subdomain_id_type, MAST::OutputFunctionBase*> output_map;
    output_map.insert(std::pair<libMesh::subdomain_id_type, MAST::OutputFunctionBase*>
                      (0, &output));
    
    
    // get reference to the element in this mesh
    const libMesh::Elem& elem = **(v._mesh->local_elements_begin());
    
    // now create the structural element
    std::unique_ptr<MAST::StructuralElementBase>
    e(MAST::build_structural_element(*v._structural_sys,
                                     elem,
                                     *v._p_card).release());
    
    
    // number of dofs in this element
    const libMesh::DofMap& dofmap = v._sys->get_dof_map();
    std::vector<unsigned int> dof_ids;
    dofmap.dof_indices(&elem, dof_ids);
    
    const unsigned int ndofs = (unsigned int)dof_ids.size();

    // make sure that the input dof vector is properly sized
    libmesh_assert_equal_to(ndofs, x0.size());
    
    // now get the residual and Jacobian evaluations
    RealVectorX
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
        
        output.clear(false);
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
            
            output.clear(false);
        }
    }
    
    // now compare the matrices
    BOOST_TEST_MESSAGE("  ** Stress Derivative  wrt State **");
    BOOST_CHECK(MAST::compare_matrix(   dstressdX_fd,  dstressdX0,    tol));
    BOOST_TEST_MESSAGE("  ** Strain Derivative  wrt State **");
    BOOST_CHECK(MAST::compare_matrix(  dstraindX_fd,   dstraindX0,    tol));
    BOOST_TEST_MESSAGE("  ** VM-Stress Derivative  wrt State **");
    BOOST_CHECK(MAST::compare_vector(     dvm_dX_fd,    dvm_dX0,      tol));
    BOOST_TEST_MESSAGE("  ** VM_func-Stress Derivative  wrt State **");
    BOOST_CHECK(MAST::compare_vector(    dvmf_dX_fd,   dvmf_dX0,      tol));
    
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
            
            output.clear(false);
        }

        
        // now the finite difference sensitivity
        // identify the perturbation in the parameter
        p0           = f();
        (fabs(p0) > 0)?  dp=delta*p0 : dp=delta;
        f()         += dp;
        
        // get the stress, strain and von Mises stress sensitivity
        e->volume_output_quantity(false, false, output_map);
        
        // reset the parameter value
        f()        = p0;
        
        // next, check the total derivative of the quantity wrt the parameter
        {
            const std::vector<MAST::StressStrainOutputBase::Data*>&
            stress_data = output.get_stress_strain_data_for_elem(&elem);
            
            libmesh_assert_equal_to(stress_data.size(), 1); // this should have one element
            
            stress              = (stress_data[0]->stress() - stress0)/dp;
            strain              = (stress_data[0]->strain() - strain0)/dp;
            vm                  = (stress_data[0]->von_Mises_stress() - vm0)/dp;
            dvmf_dp_fd          =
            (output.von_Mises_p_norm_functional_for_all_elems(pval)-vmf0)/dp;
            
            output.clear(false);
        }
        
        // compare and benchmark the values
        BOOST_TEST_MESSAGE("  ** dStress/dp (partial) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_vector(stress, dstressdp, tol));
        BOOST_TEST_MESSAGE("  ** dStrain/dp (partial) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_vector(strain, dstraindp, tol));
        BOOST_TEST_MESSAGE("  ** dVM-Stress/dp (partial) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_value (    vm,     dvmdp, tol));
        BOOST_TEST_MESSAGE("  ** dVM_func-Stress/dp (partial) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_value (dvmf_dp_fd,  dvmf_dp, tol));
    }
    
    
    //
    // now the total sensitivity including the sensitivity of the state.
    // Each dof is independently assigned a nonzero sensitivity value
    //
    for (unsigned int i=0; i<v._params_for_sensitivity.size(); i++) {
        
        MAST::Parameter& f = *v._params_for_sensitivity[i];
        
        
        // now perturb each element of the vector independently
        for (unsigned int j=0; j<ndofs; j++) {
            
            // first calculate the analytical sensitivity
            e->set_solution(x0);
            // identify the perturbation in the parameter
            // if delta=1.e-5, and dp=1.e-3, then delta/dp = 1.e-2
            // magnitude of delta/dp should not matter, since the
            // sensitivity is linear in it. We really want delta and dp to be
            // small so that finite differencing can be accurately done
            p0           = f();
            (fabs(p0) > 0)?  dp=std::max(delta*p0, delta) : dp=delta;
            // now, set the solution sensitivity for analytical sensitivity
            x            = RealVectorX::Zero(ndofs);
            x(j)         = delta/dp;
            e->set_solution(x, true);  // sets the sensitivity
            e->sensitivity_param  = &f;

            // sensitivity output
            e->volume_output_quantity(false, true, output_map);

            {
                // copy it for comparison
                const std::vector<MAST::StressStrainOutputBase::Data*>&
                stress_data = output.get_stress_strain_data_for_elem(&elem);
                
                libmesh_assert_equal_to(stress_data.size(), 1); // this should have one element
                
                dstressdp           = stress_data[0]->get_stress_sensitivity(&f);
                dstraindp           = stress_data[0]->get_strain_sensitivity(&f);
                dvmdp               = stress_data[0]->dvon_Mises_stress_dp  (&f);
                dvmf_dp             =
                output.von_Mises_p_norm_functional_sensitivity_for_all_elems(pval, &f);
                
                output.clear(false);
            }
            
            
            // now pertueb the parameter for finite differencing
            f()         += dp;
            
            // remove the sensitivity association
            e->sensitivity_param  = nullptr;

            // first perturb the solution for fd sensitivity
            x       = x0;
            x(j)   += delta;
            e->set_solution(x);
            // now, set the solution sensitivity for analytical sensitivity
            x       = RealVectorX::Zero(ndofs);
            e->set_solution(x, true);  // sets the sensitivity to zero
            
            // get the stress, strain and von Mises stress for finite-differencing
            e->volume_output_quantity(false, false, output_map);

            // reset the parameter value
            f()        = p0;

            // next, check the total derivative of the quantity wrt the parameter
            {
                const std::vector<MAST::StressStrainOutputBase::Data*>&
                stress_data = output.get_stress_strain_data_for_elem(&elem);
                
                libmesh_assert_equal_to(stress_data.size(), 1); // this should have one element
                
                stress              = (stress_data[0]->stress() - stress0)/dp;
                strain              = (stress_data[0]->strain() - strain0)/dp;
                vm                  = (stress_data[0]->von_Mises_stress() - vm0)/dp;
                dvmf_dp_fd          =
                (output.von_Mises_p_norm_functional_for_all_elems(pval)-vmf0)/dp;
                
                output.clear(false);
            }
            
            // compare and benchmark the values
            BOOST_TEST_MESSAGE("  ** dStress/dp (total: dof = " << j << ") wrt : " << f.name() << " **");
            BOOST_CHECK(MAST::compare_vector( stress, dstressdp, tol));
            BOOST_TEST_MESSAGE("  ** dStrain/dp (total: dof = " << j << ") wrt : " << f.name() << " **");
            BOOST_CHECK(MAST::compare_vector( strain, dstraindp, tol));
            BOOST_TEST_MESSAGE("  ** dVM-Stress/dp (total: dof = " << j << ") wrt : " << f.name() << " **");
            BOOST_CHECK(MAST::compare_value (  vm,   dvmdp,  tol));
            BOOST_TEST_MESSAGE("  ** dVM_func-Stress/dp (total: dof = " << j << ") wrt : " << f.name() << " **");
            BOOST_CHECK(MAST::compare_value ( dvmf_dp_fd, dvmf_dp, tol));
        }
    }
}


BOOST_AUTO_TEST_CASE   (VonMisesStress) {

    const Real
    tol      = 1.e-2;
    
    // check the accuracy of sensitivity analysis of von Mises stress,
    // and the von Mises stress functional
    // this simulates a case with 4 different stress values for an element
    
    std::unique_ptr<libMesh::Elem> elem(new libMesh::Edge2);
    MAST::Parameter f("a", 0.);
    libMesh::Point     p;
    RealVectorX
    strain = RealVectorX::Zero(6),
    stress = RealVectorX::Zero(6);

    
    Real
    stress1  = 5.e6,
    stress2  = 15.e6,
    dstress1 = -1.e8,
    dstress2 = -3.e8,
    func     = 0.,
    dfunc    = 0.,
    JxW      = 0.1;
    
    MAST::StressStrainOutputBase  output;
    
    // the four stress values
    stress(0)   =   stress1;
    output.add_stress_strain_at_qp_location(elem.get(), p, p, stress, strain, JxW);
    stress(0)   =  -stress1;
    output.add_stress_strain_at_qp_location(elem.get(), p, p, stress, strain, JxW);
    stress(0)   =   stress2;
    output.add_stress_strain_at_qp_location(elem.get(), p, p, stress, strain, JxW);
    stress(0)   =  -stress2;
    output.add_stress_strain_at_qp_location(elem.get(), p, p, stress, strain, JxW);

    // now, the stress sensitivity values
    const std::vector<MAST::StressStrainOutputBase::Data*>&
    data = output.get_stress_strain_data_for_elem(elem.get());
    
    // set the sensitivity for each stress
    stress(0)   =   dstress1;
    data[0]->set_sensitivity(&f, stress, strain);
    stress(0)   =  -dstress1;
    data[1]->set_sensitivity(&f, stress, strain);
    stress(0)   =   dstress2;
    data[2]->set_sensitivity(&f, stress, strain);
    stress(0)   =  -dstress2;
    data[3]->set_sensitivity(&f, stress, strain);
    
    // now check the vm stress value for each case
    BOOST_TEST_MESSAGE("   ** von Mises Stress ** ");
    BOOST_CHECK(MAST::compare_value(fabs(stress1),
                                    data[0]->von_Mises_stress(),
                                    tol));
    BOOST_CHECK(MAST::compare_value(fabs(stress1),
                                    data[1]->von_Mises_stress(),
                                    tol));
    BOOST_CHECK(MAST::compare_value(fabs(stress2),
                                    data[2]->von_Mises_stress(),
                                    tol));
    BOOST_CHECK(MAST::compare_value(fabs(stress2),
                                    data[3]->von_Mises_stress(),
                                    tol));
    
    BOOST_TEST_MESSAGE("   ** dvm-stress/dp **");
    BOOST_CHECK(MAST::compare_value(dstress1,
                                    data[0]->dvon_Mises_stress_dp(&f),
                                    tol));
    BOOST_CHECK(MAST::compare_value(dstress1,
                                    data[1]->dvon_Mises_stress_dp(&f),
                                    tol));
    BOOST_CHECK(MAST::compare_value(dstress2,
                                    data[2]->dvon_Mises_stress_dp(&f),
                                    tol));
    BOOST_CHECK(MAST::compare_value(dstress2,
                                    data[3]->dvon_Mises_stress_dp(&f),
                                    tol));

    BOOST_TEST_MESSAGE("   ** vm-stress functional **");
    func    =  stress2/pow(4.*JxW,0.5)*
    pow(pow(fabs( stress1)/stress2,2)*JxW +
        pow(fabs(-stress1)/stress2,2)*JxW +
        pow(fabs( stress2)/stress2,2)*JxW +
        pow(fabs(-stress2)/stress2,2)*JxW,0.5);
    
    BOOST_CHECK(MAST::compare_value(func,
                                    output.von_Mises_p_norm_functional_for_all_elems(2),
                                    tol));
    
    
    BOOST_TEST_MESSAGE("   ** dvm-stress functional/dp **");
    // update the stress values for a small perturbation
    dfunc   =  stress2/pow(4.*JxW,0.5)*0.5*
    pow(pow(fabs( stress1)/stress2,2)*JxW +
        pow(fabs(-stress1)/stress2,2)*JxW +
        pow(fabs( stress2)/stress2,2)*JxW +
        pow(fabs(-stress2)/stress2,2)*JxW,-.5) *
    2 *(fabs( stress1)/stress2*dstress1/stress2*JxW +
        fabs(-stress1)/stress2*dstress1/stress2*JxW +
        fabs( stress2)/stress2*dstress2/stress2*JxW +
        fabs(-stress2)/stress2*dstress2/stress2*JxW);

    BOOST_CHECK(MAST::compare_value
                (dfunc,
                 output.von_Mises_p_norm_functional_sensitivity_for_all_elems(2, &f),
                 tol));
    
}



BOOST_FIXTURE_TEST_SUITE  (Structural1DStressEvaluation, MAST::BuildStructural1DElem)

BOOST_AUTO_TEST_CASE   (StressLinear1DIndependentOffset) {

    this->init(false, false);
    
    RealVectorX v;

    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(1, 0, libMesh::INVALID_ELEM, v);
    check_stress(*this, v);
    
    // pure bending deformation
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(1, 1, libMesh::INVALID_ELEM, v);
    check_stress(*this, v);
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(1, 2, libMesh::INVALID_ELEM, v);
    check_stress(*this, v);
}




BOOST_AUTO_TEST_CASE   (StressNonlinear1DIndependentOffset) {
    
    this->init(false, true);
    
    RealVectorX v;
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_stress(*this, v);
}



BOOST_AUTO_TEST_CASE   (StressLinear1DDependentOffset) {
    
    this->init(true, false);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(1, 0, libMesh::INVALID_ELEM, v);
    check_stress(*this, v);
    
    // pure bending deformation
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(1, 1, libMesh::INVALID_ELEM, v);
    check_stress(*this, v);
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(1, 2, libMesh::INVALID_ELEM, v);
    check_stress(*this, v);
}



BOOST_AUTO_TEST_CASE   (StressNonlinear1DDependentOffset) {
    
    this->init(true, true);
    
    RealVectorX v;
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_stress(*this, v);
}


BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE  (Structural2DStressEvaluation, MAST::BuildStructural2DElem)

BOOST_AUTO_TEST_CASE   (StressLinear2DIndependentOffsetQUAD4) {

    this->init(false, false, libMesh::QUAD4);
    
    RealVectorX v;

    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(2, 0, libMesh::QUAD4, v);
    check_stress(*this, v);
    
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(2, 1, libMesh::QUAD4, v);
    check_stress(*this, v);
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(2, 2, libMesh::QUAD4, v);
    check_stress(*this, v);
}



BOOST_AUTO_TEST_CASE   (Stress2DNonlinearIndependentOffsetQUAD4) {
    
    this->init(false, true, libMesh::QUAD4);
    
    RealVectorX v;
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_stress(*this, v);
}



BOOST_AUTO_TEST_CASE   (StressLinear2DIndependentOffsetTRI3) {
    
    this->init(false, false, libMesh::TRI3);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(2, 0, libMesh::TRI3, v);
    check_stress(*this, v);
    
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(2, 1, libMesh::TRI3, v);
    check_stress(*this, v);
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(2, 2, libMesh::TRI3, v);
    check_stress(*this, v);
}



BOOST_AUTO_TEST_CASE   (StressNonlinear2DIndependentOffsetTRI3) {
    
    this->init(false, true, libMesh::TRI3);
    
    RealVectorX v;
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(2, 3, libMesh::TRI3, v);
    check_stress(*this, v);
}



BOOST_AUTO_TEST_CASE   (StressLinear2DDependentOffsetQUAD4) {
    
    this->init(true, false, libMesh::QUAD4);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(2, 0, libMesh::QUAD4, v);
    check_stress(*this, v);
    
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(2, 1, libMesh::QUAD4, v);
    check_stress(*this, v);
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(2, 2, libMesh::QUAD4, v);
    check_stress(*this, v);
}



BOOST_AUTO_TEST_CASE   (StressNonlinear2DDependentOffsetQUAD4) {
    
    this->init(true, true, libMesh::QUAD4);
    
    RealVectorX v;
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_stress(*this, v);
}



BOOST_AUTO_TEST_CASE   (StressLinear2DDependentOffsetTRI3) {
    
    this->init(true, false, libMesh::TRI3);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(2, 0, libMesh::TRI3, v);
    check_stress(*this, v);
    
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(2, 1, libMesh::TRI3, v);
    check_stress(*this, v);
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(2, 2, libMesh::TRI3, v);
    check_stress(*this, v);
}




BOOST_AUTO_TEST_CASE   (StressNonlinear2DDependentOffsetTRI3) {
    
    this->init(true, true, libMesh::TRI3);
    
    RealVectorX v;
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(2, 3, libMesh::TRI3, v);
    check_stress(*this, v);
}



BOOST_AUTO_TEST_SUITE_END()


