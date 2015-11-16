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


extern const Real
delta,
tol;


BOOST_FIXTURE_TEST_SUITE  (StructuralSectionPropertyEvaluation, MAST::BuildStructural1DElem)

BOOST_AUTO_TEST_CASE   (Solid1DSection) {
    

    
    
    // get reference to the element in this mesh
    const libMesh::Elem& elem = **(_mesh->local_elements_begin());
    
    // now create the structural element
    std::auto_ptr<MAST::StructuralElementBase>
    e(MAST::build_structural_element(*_structural_sys,
                                     elem,
                                     *_p_card).release());
    

    std::auto_ptr<MAST::FieldFunction<RealMatrixX > >
    mat_stiff_A  = _p_card->stiffness_A_matrix(*e);

    
    RealMatrixX
    mat0,
    dmat,
    dmatdp;
    
    Real
    p0      = 0.,
    dp      = 0.;

    
    libMesh::Point pt;
    
    // now, check the sensitivity with respect to various parameters
    
    // NOTE:
    // First, we will set the dX/dp = 0, so that the derivative of the output
    // quantity is only the partial derivative wrt the parameter.
    //

    (*mat_stiff_A)(pt, 0., mat0);

    
    for (unsigned int i=0; i<_params_for_sensitivity.size(); i++) {
        
        MAST::Parameter& f = *_params_for_sensitivity[i];

        // get the analytical sensitivity
        mat_stiff_A->derivative(MAST::PARTIAL_DERIVATIVE,
                                f,
                                pt,
                                0,
                                dmatdp);
        
        // identify the perturbation in the parameter
        p0           = f();
        (p0 > 0)?  dp=std::max(delta*p0, delta) : dp=delta;
        f()         += dp;
        
        (*mat_stiff_A)(pt, 0., dmat);
        
        // reset the parameter value
        f()        = p0;
        
        dmat    -=  mat0;
        dmat     /= dp;
        
        // compare and benchmark the values
        BOOST_CHECK(MAST::compare_matrix ( dmat, dmatdp,  tol));
    }
}

BOOST_AUTO_TEST_SUITE_END()


