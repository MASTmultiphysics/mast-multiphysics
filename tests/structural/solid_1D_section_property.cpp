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
#include "tests/base/test_comparisons.h"
#include "tests/structural/build_structural_elem_1D.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/stress_output_base.h"


// libMesh includes
#include "libmesh/dof_map.h"



template <typename ValType>
void  check_material_property (ValType& v) {
    
    const Real
    delta    = 1.e-8,
    tol      = 1.e-2;
    
    
    
    // get reference to the element in this mesh
    const libMesh::Elem& elem = **(v._mesh->local_elements_begin());
    
    // now create the structural element
    std::unique_ptr<MAST::StructuralElementBase>
    e(MAST::build_structural_element(*v._structural_sys,
                                     elem,
                                     *v._p_card).release());
    
    
    RealMatrixX
    mat0,
    dmat,
    dmatdp;
    
    Real
    p0      = 0.,
    dp      = 0.;
    
    
    libMesh::Point pt;
    
    
    
    for (unsigned int jj=0; jj<3; jj++) {
    
        
        std::unique_ptr<MAST::FieldFunction<RealMatrixX > > mat_stiff;
        
        switch (jj) {
            case 0: {
                BOOST_TEST_MESSAGE("**** Sensitivity for A matrix **");
                mat_stiff = v._p_card->stiffness_A_matrix(*e);
            }
                break;

            case 1: {
                BOOST_TEST_MESSAGE("**** Sensitivity for B matrix **");
                mat_stiff = v._p_card->stiffness_B_matrix(*e);
            }
                break;

            case 2: {
                BOOST_TEST_MESSAGE("**** Sensitivity for D matrix **");
                mat_stiff = v._p_card->stiffness_D_matrix(*e);
            }
                break;

            default:
                break;
        }
        
        // now, check the sensitivity with respect to various parameters
        
        // NOTE:
        // First, we will set the dX/dp = 0, so that the derivative of the output
        // quantity is only the partial derivative wrt the parameter.
        //
        
        (*mat_stiff)(pt, 0., mat0);
        
        
        for (unsigned int i=0; i<v._params_for_sensitivity.size(); i++) {
            
            MAST::Parameter& f = *v._params_for_sensitivity[i];
            
            // get the analytical sensitivity
            mat_stiff->derivative(  f,
                                  pt,
                                  0,
                                  dmatdp);
            
            // identify the perturbation in the parameter
            p0           = f();
            (fabs(p0) > 0)?  dp=std::max(delta*p0, delta) : dp=delta;
            f()         += dp;
            
            (*mat_stiff)(pt, 0., dmat);
            
            // reset the parameter value
            f()        = p0;
            
            dmat    -=  mat0;
            dmat     /= dp;
            
            // compare and benchmark the values
            BOOST_TEST_MESSAGE("  ** dprop/dp  wrt : " << f.name() << " **");
            BOOST_CHECK(MAST::compare_matrix ( dmat, dmatdp,  tol));
        }
    }
}


BOOST_FIXTURE_TEST_SUITE  (Structural1DSectionPropertyEvaluation,
                           MAST::BuildStructural1DElem)


BOOST_AUTO_TEST_CASE   (Property1DSensitivityIndependentOffset) {
 
    this->init(false, false); // nonlinear strain does not influence this
    check_material_property<MAST::BuildStructural1DElem>(*this);
}


BOOST_AUTO_TEST_CASE   (Property1DSensitivityDependentOffset) {
    
    this->init(true, false); // nonlinear strain does not influence this
    check_material_property<MAST::BuildStructural1DElem>(*this);
}


BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE  (Structural2DSectionPropertyEvaluation,
                           MAST::BuildStructural2DElem)


BOOST_AUTO_TEST_CASE   (Property2DSensitivityIndependentOffsetQUAD4) {
    
    this->init(false, false, libMesh::QUAD4);
    check_material_property<MAST::BuildStructural2DElem>(*this);
}


BOOST_AUTO_TEST_CASE   (Property2DSensitivityDependentOffsetQUAD4) {
    
    this->init(true, false, libMesh::QUAD4);
    check_material_property<MAST::BuildStructural2DElem>(*this);
}


BOOST_AUTO_TEST_CASE   (Property2DSensitivityIndependentOffsetTRI3) {
    
    this->init(false, false, libMesh::TRI3);
    check_material_property<MAST::BuildStructural2DElem>(*this);
}


BOOST_AUTO_TEST_CASE   (Property2DSensitivityDependentOffsetTRI3) {
    
    this->init(true, false, libMesh::TRI3);
    check_material_property<MAST::BuildStructural2DElem>(*this);
}



BOOST_AUTO_TEST_SUITE_END()


