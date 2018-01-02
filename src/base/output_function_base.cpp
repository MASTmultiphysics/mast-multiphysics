/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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

// MAST includes
#include "base/output_function_base.h"

MAST::OutputFunctionBase::OutputFunctionBase(MAST::OutputQuantityType t):
_type(t),
_eval_mode(MAST::CENTROID) {
    
}



MAST::OutputFunctionBase::~OutputFunctionBase() {
    
}



MAST::OutputQuantityType
MAST::OutputFunctionBase::type() const {
    
    return _type;
}


void
MAST::OutputFunctionBase::
set_points_for_evaluation(const std::vector<libMesh::Point>& pts) {
    
    // make sure that some points were specified
    libmesh_assert(pts.size());
    
    _eval_mode = MAST::SPECIFIED_POINTS;
    _eval_points = pts;
}



const std::vector<libMesh::Point>&
MAST::OutputFunctionBase::get_points_for_evaluation() const {

    // make sure that the data was provide
    libmesh_assert(_eval_mode == MAST::SPECIFIED_POINTS);
    return _eval_points;
}


