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


// MAST includes
#include "elasticity/stress_output_base.h"



MAST::StressStrainOutputBase::StressStrainOutputBase():
MAST::OutputFunctionBase(),
_if_evaluate_stress(true),
_if_evaluate_strain(false) {
    
}



MAST::StressStrainOutputBase::~StressStrainOutputBase() {
    
}



void
MAST::StressStrainOutputBase::set_evaluate_stress(bool val) {
    _if_evaluate_stress = val;
}



void
MAST::StressStrainOutputBase::set_evaluate_strain(bool val) {
    _if_evaluate_strain = val;
}


bool
MAST::StressStrainOutputBase::if_evaluate_stress() const {
    return _if_evaluate_stress;
}



bool
MAST::StressStrainOutputBase::if_evaluate_strain() const {
    return _if_evaluate_strain;
}



void
MAST::StressStrainOutputBase::
add_stress_at_qp_location(const libMesh::Point& quadrature_pt,
                          const libMesh::Point& physical_pt,
                          const RealVectorX& vec) {
    
}



void
MAST::StressStrainOutputBase::
add_strain_at_qp_location(const libMesh::Point& quadrature_pt,
                          const libMesh::Point& physical_pt,
                          const RealVectorX& vec) {
    
}

