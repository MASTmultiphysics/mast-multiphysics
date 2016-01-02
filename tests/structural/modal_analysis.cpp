/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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
#include "examples/structural/beam_modal_analysis/beam_modal_analysis.h"
#include "tests/base/test_comparisons.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_element_base.h"

// libMesh includes
#include "libmesh/numeric_vector.h"



BOOST_FIXTURE_TEST_SUITE  (Structural1DBeamModalAnalysis,
                           MAST::BeamModalAnalysis)

BOOST_AUTO_TEST_CASE   (BeamModalSolution) {
    
    //this->solve();
}


BOOST_AUTO_TEST_SUITE_END()



