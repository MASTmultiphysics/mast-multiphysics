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
#include <iomanip>

// MAST includes
#include "level_set/filter_base.h"

// libMesh includes
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#ifdef LIBMESH_HAVE_NANOFLANN
#include "libmesh/nanoflann.hpp"
#endif


MAST::FilterBase::FilterBase(libMesh::System& sys,
                             const Real radius,
                             const std::set<unsigned int>& dv_dof_ids):
_level_set_system  (sys),
_radius            (radius),
_level_set_fe_size (0.),
_dv_dof_ids        (dv_dof_ids) {
    
    libmesh_assert_greater(radius, 0.);
    
#ifdef LIBMESH_HAVE_NANOFLANN
    _init2();  // KD-tree search using NanoFlann
#else
    _init(); // linear filter search
#endif
}


MAST::FilterBase::~FilterBase() {
    
}


void
MAST::FilterBase::compute_filtered_values(const libMesh::NumericVector<Real>& input,
                                          libMesh::NumericVector<Real>& output,
                                          bool close_vector) const {
    
    libmesh_assert_equal_to(input.size(), _filter_map.size());
    libmesh_assert_equal_to(output.size(), _filter_map.size());
    
    output.zero();
    
    std::vector<Real> input_vals(input.size(), 0.);
    input.localize(input_vals);
    
    std::map<unsigned int, std::vector<std::pair<unsigned int, Real>>>::const_iterator
    map_it   = _filter_map.begin(),
    map_end  = _filter_map.end();
    
    for ( ; map_it != map_end; map_it++) {
        
        std::vector<std::pair<unsigned int, Real>>::const_iterator
        vec_it  = map_it->second.begin(),
        vec_end = map_it->second.end();
        
        for ( ; vec_it != vec_end; vec_it++) {
            if (map_it->first >= input.first_local_index() &&
                map_it->first <  input.last_local_index()) {
                
                if (_dv_dof_ids.count(map_it->first))
                    output.add(map_it->first, input_vals[vec_it->first] * vec_it->second);
                else
                    output.set(map_it->first, input_vals[map_it->first]);
            }
        }
    }
    
    if (close_vector)
        output.close();
}



void
MAST::FilterBase::compute_filtered_values(std::map<unsigned int, Real>& nonzero_input,
                                          libMesh::NumericVector<Real>& output,
                                          bool close_vector) const {
    
    libmesh_assert_equal_to(output.size(), _filter_map.size());
    
    output.zero();
    
    std::map<unsigned int, std::vector<std::pair<unsigned int, Real>>>::const_iterator
    map_it   = _filter_map.begin(),
    map_end  = _filter_map.end();
    
    for ( ; map_it != map_end; map_it++) {
        
        std::vector<std::pair<unsigned int, Real>>::const_iterator
        vec_it  = map_it->second.begin(),
        vec_end = map_it->second.end();
        
        for ( ; vec_it != vec_end; vec_it++) {
            if (nonzero_input.count(vec_it->first)) {
                
                if (output.type() == libMesh::SERIAL ||
                    (map_it->first >= output.first_local_index() &&
                     map_it->first <  output.last_local_index())) {
                    
                    if (_dv_dof_ids.count(map_it->first))
                        output.add(map_it->first, nonzero_input[vec_it->first] * vec_it->second);
                    else
                        output.set(map_it->first, nonzero_input[map_it->first]);
                }
            }
        }
    }
    
    if (close_vector)
        output.close();
}




void
MAST::FilterBase::compute_filtered_values(const std::vector<Real>& input,
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
        
        for ( ; vec_it != vec_end; vec_it++) {
            if (_dv_dof_ids.count(map_it->first))
                output[map_it->first] += input[vec_it->first] * vec_it->second;
            else
                output[map_it->first] += input[map_it->first];
        }
    }
}


bool
MAST::FilterBase::
if_elem_in_domain_of_influence(const libMesh::Elem& elem,
                               const libMesh::Node& level_set_node) const {
    
    Real
    d    = 1.e12; // arbitrarily large value to initialize the search
    
    libMesh::Point
    pt;
    
    // first get the smallest distance from the node to the element nodes
    for (unsigned int i=0; i<elem.n_nodes(); i++) {
        pt  = elem.point(i);
        pt -= level_set_node;
        
        if (pt.norm() < d)
            d = pt.norm();
    }
    
    // if this distance is outside the domain of influence, then this
    // element is not influenced by the design variable
    if (d > _radius + _level_set_fe_size)
        return false;
    else
        return true;
}


#ifdef LIBMESH_HAVE_NANOFLANN
// Nanoflann uses "duck typing" to allow users to define their own adaptors...
template <unsigned int Dim>
class NanoflannMeshAdaptor
{
private:
    // Constant reference to the Mesh we are adapting for use in Nanoflann
    const libMesh::MeshBase & _mesh;
    
public:
    NanoflannMeshAdaptor (const libMesh::MeshBase & mesh) :
    _mesh(mesh)
    {}
    
    /**
     * libMesh \p Point coordinate type
     */
    typedef Real coord_t;
    
    /**
     * Must return the number of data points
     */
    inline size_t
    kdtree_get_point_count() const { return _mesh.n_nodes(); }
    
    /**
     * Returns the distance between the vector "p1[0:size-1]"
     * and the data point with index "idx_p2" stored in _mesh
     */
    inline coord_t
    kdtree_distance(const coord_t * p1,
                    const size_t idx_p2,
                    size_t size) const {
        
        libmesh_assert_equal_to (size, Dim);
        
        // Construct a libmesh Point object from the input coord_t.  This
        // assumes LIBMESH_DIM==3.
        libMesh::Point point1(p1[0],
                              size > 1 ? p1[1] : 0.,
                              size > 2 ? p1[2] : 0.);
        
        // Get the referred-to point from the Mesh
        const libMesh::Point & point2 = _mesh.point(idx_p2);
        
        // Compute Euclidean distance
        return (point1 - point2).norm_sq();
    }
    
    /**
     * Returns the dim'th component of the idx'th point in the class:
     * Since this is inlined and the "dim" argument is typically an immediate value, the
     *  "if's" are actually solved at compile time.
     */
    inline coord_t
    kdtree_get_pt(const size_t idx, int dim) const
    {
        libmesh_assert_less (dim, (int) Dim);
        libmesh_assert_less (idx, _mesh.n_nodes());
        libmesh_assert_less (dim, 3);
        
        return _mesh.point(idx)(dim);
    }
    
    /**
     * Optional bounding-box computation: return false to default to a standard bbox computation loop.
     * Return true if the BBOX was already computed by the class and returned in "bb" so it can be
     * avoided to redo it again. Look at bb.size() to find out the expected dimensionality
     * (e.g. 2 or 3 for point clouds)
     */
    template <class BBOX>
    bool kdtree_get_bbox(BBOX & /* bb */) const { return false; }
};


void
MAST::FilterBase::_init2() {
    
    //libmesh_assert(!_filter_map.size());
    
    libMesh::MeshBase& mesh = _level_set_system.get_mesh();
    
    // currently implemented for replicated mesh
    libmesh_assert(mesh.is_replicated());
    
    // Loop over nodes to try and detect duplicates.  We use nanoflann
    // for this, inspired by
    // https://gist.github.com/jwpeterson/7a36f9f794df67d51126#file-detect_slit-cc-L65
    // which was inspired by nanoflann example in libMesh source:
    // contrib/nanoflann/examples/pointcloud_adaptor_example.cpp
    
    // Declare a type templated on NanoflannMeshAdaptor
    typedef nanoflann::L2_Simple_Adaptor<Real, NanoflannMeshAdaptor<3> > adatper_t;
    
    // Declare a KDTree type based on NanoflannMeshAdaptor
    typedef nanoflann::KDTreeSingleIndexAdaptor<adatper_t, NanoflannMeshAdaptor<3>, 3> kd_tree_t;
    
    // Build adaptor and tree objects
    NanoflannMeshAdaptor<3> mesh_adaptor(mesh);
    kd_tree_t kd_tree(3, mesh_adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(/*max leaf=*/10));
    
    // Construct the tree
    kd_tree.buildIndex();
    
    Real
    d_12 = 0.,
    sum  = 0.;
    
    unsigned int
    dof_1,
    dof_2;
    
    libMesh::MeshBase::const_node_iterator
    node_it      =  mesh.nodes_begin(),
    node_end     =  mesh.nodes_end();
    
    // For every node in the mesh, search the KDtree and find any
    // nodes at _radius distance from the current
    // node being searched... this will be added to the .
    for (; node_it != node_end; node_it++) {
        
        const libMesh::Node* node = *node_it;
        
        dof_1 = node->dof_number(_level_set_system.number(), 0, 0);
        
        Real query_pt[3] = {(*node)(0), (*node)(1), (*node)(2)};
        
        std::vector<std::pair<size_t, Real>>
        indices_dists;
        nanoflann::RadiusResultSet<Real, size_t>
        resultSet(_radius*_radius, indices_dists);
        
        kd_tree.findNeighbors(resultSet, query_pt, nanoflann::SearchParams());
        
        sum       = 0.;
        
        for (unsigned r=0; r<indices_dists.size(); ++r) {
            
            d_12 = std::sqrt(indices_dists[r].second);
            
            // the distance of this node should be less than or equal to the
            // specified search radius
            libmesh_assert_less_equal(d_12, _radius);
            
            sum  += _radius - d_12;
            dof_2 = mesh.node_ptr(indices_dists[r].first)->dof_number(_level_set_system.number(), 0, 0);
            
            _filter_map[dof_1].push_back(std::pair<unsigned int, Real>(dof_2, _radius - d_12));
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
    
    // compute the largest element size
    libMesh::MeshBase::const_element_iterator
    e_it          = mesh.elements_begin(),
    e_end         = mesh.elements_end();
    
    for ( ; e_it != e_end; e_it++) {
        const libMesh::Elem* e = *e_it;
        d_12 = e->hmax();
        
        if (_level_set_fe_size < d_12)
            _level_set_fe_size = d_12;
    }
}
#endif


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
    
    // compute the largest element size
    libMesh::MeshBase::const_element_iterator
    e_it          = mesh.elements_begin(),
    e_end         = mesh.elements_end();
    
    for ( ; e_it != e_end; e_it++) {
        const libMesh::Elem* e = *e_it;
        d_12 = e->hmax();
        
        if (_level_set_fe_size < d_12)
            _level_set_fe_size = d_12;
    }
}


void
MAST::FilterBase::print(std::ostream& o) const {
    
    o << "Filter radius: " << _radius << std::endl;
    
    o
    << std::setw(20) << "Filtered ID"
    << std::setw(20) << "Dependent Vars" << std::endl;
    
    std::map<unsigned int, std::vector<std::pair<unsigned int, Real>>>::const_iterator
    map_it   = _filter_map.begin(),
    map_end  = _filter_map.end();
    
    for ( ; map_it != map_end; map_it++) {
        
        o
        << std::setw(20) << map_it->first;
        
        std::vector<std::pair<unsigned int, Real>>::const_iterator
        vec_it  = map_it->second.begin(),
        vec_end = map_it->second.end();
        
        for ( ; vec_it != vec_end; vec_it++) {
            
            if (_dv_dof_ids.count(map_it->first))
                o
                << " : " << std::setw(8) << vec_it->first
                << " (" << std::setw(8) << vec_it->second << " )";
            else
                libMesh::out << " : " << map_it->first;
        }
        libMesh::out << std::endl;
    }
}

