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
#include "level_set/sub_elem_mesh_refinement.h"
#include "level_set/level_set_intersection.h"
#include "level_set/sub_elem_node_map.h"

// libMesh includes
#include "libmesh/mesh_communication.h"
#include "libmesh/partitioner.h"
#include "libmesh/dof_map.h"
#include "libmesh/boundary_info.h"


MAST::SubElemMeshRefinement::SubElemMeshRefinement(libMesh::MeshBase& mesh,
                                                   libMesh::System&   sys):
libMesh::System::Constraint          (),
_initialized                         (false),
_strong_discontinuity                (false),
_negative_level_set_subdomain_offset (0),
_inactive_subdomain_offset           (0),
_level_set_boundary_id               (0),
_mesh                                (mesh),
_system                              (sys),
_node_map                            (new MAST::SubElemNodeMap) {
    
}



MAST::SubElemMeshRefinement::~SubElemMeshRefinement() {

    this->clear_mesh();
    delete _node_map;
}



bool
MAST::SubElemMeshRefinement::process_mesh(const MAST::FieldFunction<Real>& phi,
                                          bool strong_discontinuity,
                                          Real time,
                                          unsigned int negative_level_set_subdomain_offset,
                                          unsigned int inactive_subdomain_offset,
                                          unsigned int level_set_boundary_id) {

    libmesh_assert(!_initialized);
    
    // currently only implemented for replicated mesh
    libmesh_assert(_mesh.is_replicated());
    
    // if strong discontinuity is required then coincident nodes are
    // created which will need unique_id support in libMesh
    if (strong_discontinuity) {
#ifndef LIBMESH_ENABLE_UNIQUE_ID
        libmesh_assert(false);
#endif
    }
    
    MAST::LevelSetIntersection intersect;
    
    libMesh::MeshBase::element_iterator
    e_it    =  _mesh.active_local_elements_begin(),
    e_end   =  _mesh.active_local_elements_end();

    // for a replicated mesh all processors have to do the exct same operations
    // to the mesh in the same order. Hence, we modify the iterators.
    if (_mesh.is_replicated()) {
        
        e_it    =  _mesh.active_elements_begin(),
        e_end   =  _mesh.active_elements_end();
    }

    // first we need to identify all the elements that will be refined.
    // Then we will iterate over all of them. Otherwise, the addition of
    // new elemnts can invalidate the element iterators.
    std::vector<libMesh::Elem*>
    elems_to_partition;
    
    for ( ; e_it != e_end; e_it++) {
        
        libMesh::Elem* elem = *e_it;
        
        intersect.init(phi, *elem, time, _mesh.max_elem_id(), _mesh.max_node_id());

        if (intersect.if_intersection_through_elem() ||
            ((intersect.get_intersection_mode() == MAST::COLINEAR_EDGE ||
              intersect.get_intersection_mode() == MAST::THROUGH_NODE) &&
             intersect.get_sub_elems_negative_phi().size() == 1) ||
            intersect.if_elem_on_negative_phi())
            elems_to_partition.push_back(elem);
        
        intersect.clear();
    }
    
    
    // now we process only the selected elements
    bool
    mesh_changed = false;
    
    for (unsigned int i=0; i<elems_to_partition.size(); i++) {
        
        libMesh::Elem* elem = elems_to_partition[i];
        
        intersect.init(phi, *elem, time, _mesh.max_elem_id(), _mesh.max_node_id());
        
        // if the intersection is through the element then
        if (intersect.if_intersection_through_elem()) {
            
            _process_sub_elements(strong_discontinuity,
                                  negative_level_set_subdomain_offset,
                                  level_set_boundary_id,
                                  *elem,
                                  intersect,
                                  true,
                                  intersect.get_sub_elems_positive_phi());
            
            _process_sub_elements(strong_discontinuity,
                                  negative_level_set_subdomain_offset,
                                  level_set_boundary_id,
                                  *elem,
                                  intersect,
                                  false,
                                  intersect.get_sub_elems_negative_phi());
            
            // since the element has been partitioned, we set its subdomain
            // id with an offset so that the assembly routines can choose to
            // ignore them
            _old_elems.push_back(std::make_pair(elem, elem->subdomain_id()));
            elem->subdomain_id() += inactive_subdomain_offset;
            mesh_changed = true;
        }
        else if ((intersect.get_intersection_mode() == MAST::COLINEAR_EDGE ||
                  intersect.get_intersection_mode() == MAST::THROUGH_NODE) &&
                 intersect.get_sub_elems_negative_phi().size() == 1) {

            if (strong_discontinuity) {
                _process_negative_element(negative_level_set_subdomain_offset,
                                          level_set_boundary_id,
                                          *elem,
                                          intersect);
                
                _old_elems.push_back(std::make_pair(elem, elem->subdomain_id()));
                elem->subdomain_id() += inactive_subdomain_offset;
            }
            else {
                
                _old_elems.push_back(std::make_pair(elem, elem->subdomain_id()));
                elem->subdomain_id() += negative_level_set_subdomain_offset;

            }
        }
        else if (intersect.if_elem_on_negative_phi()) {
            // if the element has no positive region, then we set its
            // subdomain id to that of negative level set offset
            
            _old_elems.push_back(std::make_pair(elem, elem->subdomain_id()));
            elem->subdomain_id() += negative_level_set_subdomain_offset;
        }
        
        intersect.clear();
    }

    // If the mesh changed on any processor, it changed globally
    _mesh.comm().max(mesh_changed);
    
    // And we may need to update DistributedMesh values reflecting the changes
    if (mesh_changed)
        _mesh.update_parallel_id_counts();
    
    if (mesh_changed && !_mesh.is_replicated())
    {
        libMesh::MeshCommunication().make_elems_parallel_consistent (_mesh);
        libMesh::MeshCommunication().make_new_nodes_parallel_consistent (_mesh);
#ifdef DEBUG
        _mesh.libmesh_assert_valid_parallel_ids();
#endif
    }
    
    // If we're refining a Replicated_mesh, then we haven't yet assigned
    // node processor ids.  But if we're refining a partitioned
    // Replicated_mesh, then we *need* to assign node processor ids.
    if (mesh_changed && _mesh.is_replicated() &&
        (_mesh.unpartitioned_elements_begin() ==
         _mesh.unpartitioned_elements_end()))
        libMesh::Partitioner::set_node_processor_ids(_mesh);
    
    if (mesh_changed)
        _mesh.prepare_for_use(/*skip_renumber =*/ false);
    
    _strong_discontinuity                = strong_discontinuity;
    _negative_level_set_subdomain_offset = negative_level_set_subdomain_offset;
    _inactive_subdomain_offset           = inactive_subdomain_offset;
    _level_set_boundary_id               = level_set_boundary_id;
    _initialized                         = true;
    
    return mesh_changed;
}


bool
MAST::SubElemMeshRefinement::clear_mesh() {
    
    // clear the data structure
    _hanging_node.clear();
    
    // modify the original element subdomain
    for (unsigned int i=0; i<_old_elems.size(); i++)
        _old_elems[i].first->subdomain_id() = _old_elems[i].second;
    _old_elems.clear();
    
    // remove all the new nodes and elements from the mesh
    for (unsigned int i=0; i<_new_elems.size(); i++) {
        
        // Remove this element from any neighbor
        // lists that point to it.
        _new_elems[i]->nullify_neighbors();
        
        // Remove any boundary information associated
        // with this element
        _mesh.get_boundary_info().remove(_new_elems[i]);

        _mesh.delete_elem(_new_elems[i]);
    }
    
    for (unsigned int i=0; i<_new_nodes.size(); i++)
        _mesh.delete_node(_new_nodes[i]);
    
    
    bool
    mesh_changed = false;
    
    if (_new_elems.size() || _new_nodes.size()) mesh_changed = true;
    
    _mesh.comm().max(mesh_changed);
    
    if (mesh_changed) {
        
        _mesh.update_parallel_id_counts();
        _new_elems.clear();
        _new_nodes.clear();
        _node_map->clear();
        _mesh.prepare_for_use();
    }

    if (mesh_changed && !_mesh.is_serial()) {
        
        libMesh::MeshCommunication().make_nodes_parallel_consistent (_mesh);
        
#ifdef DEBUG
        MeshTools::libmesh_assert_valid_procids<Node>(_mesh);
#endif
    }

    _negative_level_set_ids.clear();
    _strong_discontinuity                = false;
    _negative_level_set_subdomain_offset = 0;
    _inactive_subdomain_offset           = 0;
    _level_set_boundary_id               = 0;
    _initialized                         = false;
    
    return mesh_changed;
}



void
MAST::SubElemMeshRefinement::constrain() {
    
    // we will constrain only if the mesh has been processed for the level set
    // no constraint is to be added for weak discontinuity.
    if (!_initialized ||  !_strong_discontinuity) return;

    // For strong discontinuity, iterate over all elements and constrain all
    // dofs on them
    libMesh::DofMap& dof_map = _system.get_dof_map();
    
    libMesh::MeshBase::const_element_iterator       el     =
    _system.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    _system.get_mesh().active_local_elements_end();
    
    const libMesh::dof_id_type
    first_dof  = dof_map.first_dof(_mesh.comm().rank()),
    last_dof   = dof_map.end_dof(_mesh.comm().rank());

    std::vector<libMesh::dof_id_type>
    dof_indices;

    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;

        if (_negative_level_set_ids.count(elem->subdomain_id())) {
            
            dof_indices.clear();
            dof_map.dof_indices(elem, dof_indices);
            
            // constrain all dofs if they have not already been constrained.
            for (unsigned int i=0; i<dof_indices.size(); i++) {
                
                if ((dof_indices[i] >=   first_dof  || dof_indices[i] <  last_dof) &&
                    !dof_map.is_constrained_dof(dof_indices[i])) {
                    
                    libMesh::DofConstraintRow c_row;
                    dof_map.add_constraint_row(dof_indices[i], c_row, true);
                }
            }
        }
    }
    
    // now add constriant for the nodes identified as hanging node
    std::set<std::pair<const libMesh::Node*, std::pair<const libMesh::Node*, const libMesh::Node*>>>::iterator
    n_it  = _hanging_node.begin(),
    n_end = _hanging_node.end();
    
    libMesh::Point
    p;
    
    Real
    d = 0.;
    
    unsigned int
    dof_b_node1,
    dof_b_node2,
    dof_node;
    
    for ( ; n_it != n_end; n_it++) {
        
        const std::pair<const libMesh::Node*, std::pair<const libMesh::Node*, const libMesh::Node*>>
        &v   =  *n_it;

        // obtain the fraction of the node from the bounding nodes
        // distance of node from first
        p  = *v.first - *v.second.first;
        d  = p.norm();
        // distance between the bounding nodes
        p  = *v.second.second - *v.second.first;
        d  /= p.norm();
        

        // we iterate over the variables in the system
        // and add a constraint for each node
        for (unsigned int i=0; i<_system.n_vars(); i++) {
            
            // identify the dofs for the ith variable on each node
            // first, the node for which the constraint will be added
            dof_indices.clear();
            dof_map.dof_indices(v.first, dof_indices, i);
            libmesh_assert_equal_to(dof_indices.size(), 1);
            dof_node = dof_indices[0];

            // next, the first bounding node
            dof_indices.clear();
            dof_map.dof_indices(v.second.first, dof_indices, i);
            libmesh_assert_equal_to(dof_indices.size(), 1);
            dof_b_node1 = dof_indices[0];

            // next, the second bounding node
            dof_indices.clear();
            dof_map.dof_indices(v.second.second, dof_indices, i);
            libmesh_assert_equal_to(dof_indices.size(), 1);
            dof_b_node2 = dof_indices[0];

            // now create and add the constraint
            if ((dof_node >=   first_dof  || dof_node <  last_dof) &&
                !dof_map.is_constrained_dof(dof_node)) {
                
                // the constraint assumes linear variation of the value
                // between the bounding nodes
                // the constraint reads
                //      un = (1-d) ub1 + d ub2
                // or,  un - (1-d) ub1 - d ub2 = 0
                libMesh::DofConstraintRow c_row;
                c_row[dof_b_node1] = (1.-d);
                c_row[dof_b_node2] =     d;
                
                dof_map.add_constraint_row(dof_node, c_row, true);
            }
        }
    }
}



void
MAST::SubElemMeshRefinement::_process_sub_elements(bool strong_discontinuity,
                                                   unsigned int negative_level_set_subdomain_offset,
                                                   unsigned level_set_boundary_id,
                                                   libMesh::Elem& e,
                                                   MAST::LevelSetIntersection& intersect,
                                                   bool positive_phi,
                                                   const std::vector<const libMesh::Elem*>& elems) {

    libmesh_assert(!_initialized);
    
    std::map<const libMesh::Node*, const libMesh::Node*>
    intersection_object_to_mesh_node_map;

    std::vector<std::tuple<const libMesh::Node*, libMesh::Elem*, unsigned int>>
    interior_node_association;
    
    for (unsigned int i=0; i<elems.size(); i++) {
        
        const libMesh::Elem
        *sub_e = elems[i];
        
        // in case of intersection through node or collinear edge the
        // parent element is the subelement. In case of a weak discontinuity
        // nothing needs to be done since the C0 continuity will be maintained
        // due to same nodes used across positive and negative level set values.
        // However, for a strong discontinuity, the negative and positive
        // values of level set reside on separate nodes. So, for intersection
        // through a node we need to create a new node at this location.
        // Similarly, for colinear edge we need to create new nodes at
        // all nodes along this edge.
        //
        // For such cases, we will not do anythign for elements on positive
        //  level set value. Instead, we will add new elements for elements
        //  on the negative level set
        if ((sub_e == &e) && positive_phi) continue;
        
        libMesh::Elem
        *child = libMesh::Elem::build(sub_e->type()).release();
        
        // set nodes for this child element
        for (unsigned int j=0; j<sub_e->n_nodes(); j++) {
            
            const libMesh::Node
            *sub_e_node = sub_e->node_ptr(j);
            
            if (!intersect.if_node_is_new(*sub_e_node)) {
                
                // this is a node from the parent element. So, we use is as is
                child->set_node(j) = const_cast<libMesh::Node*>(sub_e_node);
                
                // keep track of nodes for the addition of interior nodes
                intersection_object_to_mesh_node_map[sub_e_node] = sub_e_node;
            }
            else if (intersect.if_interior_node(*sub_e_node)) {
                
                // this node will be added after all edge nodes have been
                // added. This will be added as the jth node of the
                // child elem
                interior_node_association.push_back
                (std::tuple<const libMesh::Node*, libMesh::Elem*, unsigned int>(sub_e_node, child, j));
            }
            else {
                
                // this is a new node. So, we ask the intersection object
                // about the bounding nodes and add them to this new node
                std::pair<const libMesh::Node*, const libMesh::Node*>
                bounding_nodes = intersect.get_bounding_nodes_for_node(*sub_e_node);
                
                libMesh::Node*
                child_node = _add_node(*sub_e_node,
                                       strong_discontinuity,
                                       positive_phi,
                                       e.processor_id(),
                                       bounding_nodes);
                child->set_node(j) = child_node;
                
                // identify this node to as a hanging node
                if (intersect.if_hanging_node(sub_e_node))
                    _hanging_node.insert(std::make_pair(child_node, bounding_nodes));
                    
                // keep track for nodes for the addition of interior nodes
                intersection_object_to_mesh_node_map[sub_e_node] = child_node;
            }
        }
        
        // set flags for this child element
        //child->set_refinement_flag(libMesh::Elem::JUST_REFINED);
        child->set_p_level(e.p_level());
        child->set_p_refinement_flag(e.p_refinement_flag());
        child->set_n_systems(e.n_systems());
        // we need to offset the subdomain id for an element type that is
        // not the same as the parent element since exodus output requires
        // different subdomain ids for different element types.
        if (child->type() == e.type())
            child->subdomain_id() = e.subdomain_id();
        else
            child->subdomain_id() = e.subdomain_id()+1;
        
        // the negative level set is offset by the specified value so that
        // the assembly routines can deal with them separately
        if (!positive_phi) {
            
            child->subdomain_id() += negative_level_set_subdomain_offset;
            _negative_level_set_ids.insert(child->subdomain_id());
        }
        
#if (LIBMESH_MAJOR_VERSION == 1 && LIBMESH_MINOR_VERSION >= 5)
        libmesh_assert_equal_to (child->n_extra_integers(),
                                 e.n_extra_integers());
        for (unsigned int j=0; j != e.n_extra_integers(); ++j)
            child->set_extra_integer(j, e.get_extra_integer(j));
#endif
        _mesh.add_elem(child);

        if (intersect.has_side_on_interface(*sub_e))
            _mesh.boundary_info->add_side(child,
                                          intersect.get_side_on_interface(*sub_e),
                                          level_set_boundary_id);
        
        _new_elems.push_back(child);
    }
    
    
    // now process the interior nodes
    for (unsigned int i=0; i<interior_node_association.size(); i++) {

        const libMesh::Node
        *sub_e_node = std::get<0>(interior_node_association[i]);

        libMesh::Elem
        *child = std::get<1>(interior_node_association[i]);

        unsigned int
        node_num = std::get<2>(interior_node_association[i]);
        
        std::pair<const libMesh::Node*, const libMesh::Node*>
        bounding_nodes = intersect.get_bounding_nodes_for_node(*sub_e_node);
        
        // since the nodes in the map are based on newly created nodes, and
        // not those created by the intersection object, we replace the
        // bounding_nodes pair with those that were created here.
        libmesh_assert(intersection_object_to_mesh_node_map.count(bounding_nodes.first));
        libmesh_assert(intersection_object_to_mesh_node_map.count(bounding_nodes.second));
        bounding_nodes.first  = intersection_object_to_mesh_node_map[bounding_nodes.first];
        bounding_nodes.second = intersection_object_to_mesh_node_map[bounding_nodes.second];
        
        libMesh::Node*
        child_node = _add_node(*sub_e_node,
                               strong_discontinuity,
                               positive_phi,
                               e.processor_id(),
                               bounding_nodes);
        
        child->set_node(node_num) = child_node;
    }
}




void
MAST::SubElemMeshRefinement::_process_negative_element(unsigned int negative_level_set_subdomain_offset,
                                                       unsigned level_set_boundary_id,
                                                       libMesh::Elem& e,
                                                       MAST::LevelSetIntersection& intersect) {
    
    libmesh_assert(!_initialized);
    
    std::set<libMesh::Node*>
    side_nodes;
    
    // get the nodes on the side with the interface that will be replaced
    // with new nodes on the negative phi
    if (intersect.get_intersection_mode() == MAST::THROUGH_NODE)
        side_nodes.insert(e.node_ptr(intersect.node_on_boundary()));
    else if (intersect.get_intersection_mode() == MAST::COLINEAR_EDGE) {
        
        unsigned int i = intersect.edge_on_boundary();
        std::unique_ptr<libMesh::Elem> side(e.side_ptr(i));
        for (unsigned int j=0; j<side->n_nodes(); j++)
            side_nodes.insert(side->node_ptr(j));
    }
    else {
        // should not get here
        libmesh_assert(false);
    }
    
    
    libMesh::Elem
    *child = libMesh::Elem::build(e.type()).release();
    
    // set nodes for this new element
    for (unsigned int j=0; j<e.n_nodes(); j++) {
        
        libMesh::Node
        *e_node = e.node_ptr(j);
        
        if (!side_nodes.count(e_node))
            child->set_node(j) = e_node;
        else {

            std::pair<const libMesh::Node*, const libMesh::Node*>
            bounding_nodes = std::make_pair(e_node, e_node);
            
            libMesh::Node*
            child_node = _add_node(*e_node,
                                   true,    // this method only deals with strong discontinuity
                                   false,   // and with nodes on negative level set
                                   e.processor_id(),
                                   bounding_nodes);
            child->set_node(j) = child_node;
        }
    }
    
    // set flags for this child element
    //child->set_refinement_flag(libMesh::Elem::JUST_REFINED);
    child->set_p_level(e.p_level());
    child->set_p_refinement_flag(e.p_refinement_flag());
    child->set_n_systems(e.n_systems());
    // we need to offset the subdomain id for an element type that is
    // not the same as the parent element since exodus output requires
    // different subdomain ids for different element types.
    if (child->type() == e.type())
        child->subdomain_id() = e.subdomain_id() + negative_level_set_subdomain_offset;
    else
        child->subdomain_id() = e.subdomain_id()+1 + negative_level_set_subdomain_offset;
    
    _negative_level_set_ids.insert(child->subdomain_id());
    
#if (LIBMESH_MAJOR_VERSION == 1 && LIBMESH_MINOR_VERSION >= 5)
    libmesh_assert_equal_to (child->n_extra_integers(),
                             e.n_extra_integers());

    for (unsigned int j=0; j != e.n_extra_integers(); ++j)
        child->set_extra_integer(j, e.get_extra_integer(j));
#endif
    _mesh.add_elem(child);
    
    _new_elems.push_back(child);
}



libMesh::Node*
MAST::SubElemMeshRefinement::_add_node(const libMesh::Point& p,
                                       bool strong_disontinuity,
                                       bool positive_phi,
                                       unsigned int processor_id,
                                       const std::pair<const libMesh::Node*, const libMesh::Node*>& bounding_nodes) {
    
    libmesh_assert(!_initialized);
    
    unsigned int
    id1 = std::min(bounding_nodes.first->id(), bounding_nodes.second->id()),
    id2 = std::max(bounding_nodes.first->id(), bounding_nodes.second->id());
    
    std::pair<libMesh::Node*, libMesh::Node*>&
    node_pair = _node_map->add(id1, id2);

    // if a weak discontinuity is requested, and if the node has already been
    // created, then make sure that both nodes are the same
    if (!strong_disontinuity) {
        
        if (node_pair.first) {
            libmesh_assert_equal_to(node_pair.first, node_pair.second);
            return node_pair.first;
        }
        else {

            // for a weak discontinuity nodes on either side of the discontinuity
            // are the same
            node_pair.first = _mesh.add_point(p, libMesh::DofObject::invalid_id, processor_id);
            _new_nodes.push_back(node_pair.first);
            
            libmesh_assert(node_pair.first);

            node_pair.first->processor_id() = libMesh::DofObject::invalid_processor_id;
            node_pair.second = node_pair.first;
            
            return node_pair.first;
        }
    }
    else {
        
        if (positive_phi) {
            
            if (node_pair.first)
                return node_pair.first;
            else {
                
                // create and store a separate node for the positive side of the
                // level set
                node_pair.first = _mesh.add_point(p, libMesh::DofObject::invalid_id, processor_id);
                _new_nodes.push_back(node_pair.first);

                libmesh_assert(node_pair.first);
                
                node_pair.first->processor_id() = libMesh::DofObject::invalid_processor_id;
                
                return node_pair.first;
            }
        }
        else {         // negative phi

            if (node_pair.second)
                return node_pair.second;
            else {
                
                // create and store a separate node for the positive side of the
                // level set
                node_pair.second = _mesh.add_point(p, libMesh::DofObject::invalid_id, processor_id);
                _new_nodes.push_back(node_pair.second);

                libmesh_assert(node_pair.second);
                
                node_pair.second->processor_id() = libMesh::DofObject::invalid_processor_id;
                
                return node_pair.second;
            }
        }
    }
}


