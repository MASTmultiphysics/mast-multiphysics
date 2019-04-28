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

// MAST includes
#include "level_set/filter_base.h"

// libMesh includes
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"



MAST::FilterBase::FilterBase(libMesh::System& sys,
                             const Real radius):
_level_set_system  (sys),
_radius            (radius) {
    
    libmesh_assert_greater(radius, 0.);
    
    _init();
}


MAST::FilterBase::~FilterBase() {
    
}


#include <iomanip>

template <>
void
MAST::FilterBase::
compute_filtered_values<libMesh::NumericVector<Real>>(const libMesh::NumericVector<Real>& input,
                                                      libMesh::NumericVector<Real>& output) const {
    
    libmesh_assert_equal_to(input.size(), _filter_map.size());
    libmesh_assert_equal_to(output.size(), _filter_map.size());
    
    output.zero();
    
    std::map<unsigned int, std::vector<std::pair<unsigned int, Real>>>::const_iterator
    map_it   = _filter_map.begin(),
    map_end  = _filter_map.end();
    
    for ( ; map_it != map_end; map_it++) {
        
        std::vector<std::pair<unsigned int, Real>>::const_iterator
        vec_it  = map_it->second.begin(),
        vec_end = map_it->second.end();
        
        for ( ; vec_it != vec_end; vec_it++)
            output.add(map_it->first, input.el(vec_it->first) * vec_it->second);
    }
    
    output.close();
}




template <>
void
MAST::FilterBase::
compute_filtered_values<std::vector<Real>>(const std::vector<Real>& input,
                                           std::vector<Real>& output) const {
    
    libmesh_assert_equal_to(input.size(), _filter_map.size());
    libmesh_assert_equal_to(output.size(), _filter_map.size());

    std::fill(output.begin(), output.end(), 0.);
    
    std::map<unsigned int, std::vector<std::pair<unsigned int, Real>>>::const_iterator
    map_it   = _filter_map.begin(),
    map_end  = _filter_map.end();
    
    for ( ; map_it != map_end; map_it++) {
        
        std::vector<std::pair<unsigned int, Real>>::const_iterator
        vec_it  = map_it->second.begin(),
        vec_end = map_it->second.end();
        
        for ( ; vec_it != vec_end; vec_it++)
            output[map_it->first] += input[vec_it->first] * vec_it->second;
    }
}



void
MAST::FilterBase::_init() {
    
    libmesh_assert(!_filter_map.size());
    
    libMesh::MeshBase& mesh = _level_set_system.get_mesh();
    
    // currently implemented for replicated mesh
    libmesh_assert(mesh.is_replicated());
    
    // iterate over all nodes to compute the
    libMesh::MeshBase::const_node_iterator
    node_it_1    =  mesh.nodes_begin(),
    node_it_2    =  mesh.nodes_begin(),
    node_end     =  mesh.nodes_end();
    
    libMesh::Point
    d;
    
    Real
    d_12 = 0.,
    sum  = 0.;
    
    unsigned int
    dof_1,
    dof_2;
    
    for ( ; node_it_1 != node_end; node_it_1++) {
        
        dof_1 = (*node_it_1)->dof_number(_level_set_system.number(), 0, 0);
        
        node_it_2 = mesh.nodes_begin();
        sum       = 0.;
        
        for ( ; node_it_2 != node_end; node_it_2++) {
            
            // compute the distance between the two nodes
            d    = (**node_it_1) - (**node_it_2);
            d_12 = d.norm();
            
            // if the nodes is within the filter radius, add it to the map
            if (d_12 <= _radius) {
                
                sum  += _radius - d_12;
                dof_2 = (*node_it_2)->dof_number(_level_set_system.number(), 0, 0);

                _filter_map[dof_1].push_back(std::pair<unsigned int, Real>(dof_2, _radius - d_12));
            }
        }
        
        libmesh_assert_greater(sum, 0.);
        
        // with the coefficients computed for dof_1, divide each coefficient
        // with the sum
        std::vector<std::pair<unsigned int, Real>>& vec = _filter_map[dof_1];
        for (unsigned int i=0; i<vec.size(); i++) {
            
            vec[i].second /= sum;
            libmesh_assert_less_equal(vec[i].second, 1.);
        }
    }
}
