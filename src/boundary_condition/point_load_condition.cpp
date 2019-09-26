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


// MAST includes
#include "boundary_condition/point_load_condition.h"


MAST::PointLoadCondition::
PointLoadCondition(MAST::BoundaryConditionType t):
MAST::BoundaryConditionBase(t) {
    
    // make sure that the type is one of the following two types
    libmesh_assert(t == MAST::POINT_LOAD ||
                   t == MAST::POINT_MOMENT);
}



MAST::PointLoadCondition::~PointLoadCondition() { }



void
MAST::PointLoadCondition::add_node(const libMesh::Node& nd) {
    
    _nodes.insert(&nd);
}


const std::set<const libMesh::Node*>&
MAST::PointLoadCondition::get_nodes() const {
    
    return _nodes;
}


std::set<const libMesh::Node*>&
MAST::PointLoadCondition::get_nodes() {
    
    return _nodes;
}



MAST::PointLoad::PointLoad(MAST::Parameter magnitude, RealVectorX direction):
MAST::FieldFunction<RealVectorX>("load"), _magnitude(magnitude),
_direction(direction) {}


MAST::PointLoad::~PointLoad() { }

void MAST::PointLoad::operator()(const libMesh::Point& p, const Real t, 
                                 RealVectorX& v) const
{
    v = _magnitude() * _direction;
}


void MAST::PointLoad::derivative(const MAST::FunctionBase& f, 
                                 const libMesh::Point& p, const Real t, 
                                 RealVectorX& v) const 
{
    v = RealVectorX::Zero(6);
    if (&f == &_magnitude)
    {
        for (uint64_t i=0; i<_direction.size(); i++)
        {
            v(i) = _direction(i);
        }
    }
}
