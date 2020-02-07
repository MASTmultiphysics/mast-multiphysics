/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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

// C++ includes
#include <vector>

// MAST includes
#include "elasticity/bending_operator.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/bernoulli_bending_operator.h"
#include "elasticity/dkt_bending_operator.h"
#include "elasticity/mindlin_bending_operator.h"
#include "elasticity/timoshenko_bending_operator.h"



MAST::BendingOperator::BendingOperator(MAST::StructuralElementBase& elem):
_structural_elem(elem),
_elem(_structural_elem.elem())
{ }



MAST::BendingOperator::~BendingOperator()
{ }



std::unique_ptr<MAST::BendingOperator1D>
MAST::build_bending_operator_1D(MAST::BendingOperatorType type,
                                MAST::StructuralElementBase& elem,
                                const std::vector<libMesh::Point>& pts) {
    
    std::unique_ptr<MAST::BendingOperator1D> rval;
    
    switch (type) {
        case MAST::BERNOULLI:
            rval.reset(new MAST::BernoulliBendingOperator(elem, pts));
            break;
            
        case MAST::TIMOSHENKO:
            rval.reset(new MAST::TimoshenkoBendingOperator(elem));
            break;
            
        case MAST::NO_BENDING:
            // nothing to be done
            break;
            
        default:
            libmesh_error(); // should not get here
            break;
    }
    
    return rval;
}


std::unique_ptr<MAST::BendingOperator2D>
MAST::build_bending_operator_2D(MAST::BendingOperatorType type,
                                MAST::StructuralElementBase& elem,
                                const std::vector<libMesh::Point>& pts) {
    
    std::unique_ptr<MAST::BendingOperator2D> rval;
    
    switch (type) {

        case MAST::DKT:
            rval.reset(new MAST::DKTBendingOperator(elem, pts));
            break;
            
        case MAST::MINDLIN:
            rval.reset(new MAST::MindlinBendingOperator(elem));
            break;
            
        case MAST::NO_BENDING:
            // nothing to be done
            break;
            
        default:
            libmesh_error(); // should not get here
            break;
    }
    
    return rval;
}
