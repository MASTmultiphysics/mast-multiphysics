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
#include "base/output_assembly_elem_operations.h"


MAST::OutputAssemblyElemOperations::OutputAssemblyElemOperations():
MAST::AssemblyElemOperations(),
_if_evaluate_on_all_elems(false) {
    
}


MAST::OutputAssemblyElemOperations::~OutputAssemblyElemOperations() {
    
}


    

void
MAST::OutputAssemblyElemOperations::
set_participating_subdomains(const std::set<libMesh::subdomain_id_type>& sids) {
    
    // this object should be in a clean state before this method is called
    libmesh_assert(!_sub_domain_ids.size());
    libmesh_assert(!_elem_subset.size());
    
    _sub_domain_ids           = sids;
    _if_evaluate_on_all_elems = false;
}
    



void
MAST::OutputAssemblyElemOperations::set_participating_elements_to_all() {
    
    // this object should be in a clean state before this method is called
    libmesh_assert(!_sub_domain_ids.size());
    libmesh_assert(!_elem_subset.size());
    
    _if_evaluate_on_all_elems = true;
}


void
MAST::OutputAssemblyElemOperations::
set_participating_elements(const std::set<const libMesh::Elem*>& elems) {
    
    // this object should be in a clean state before this method is called
    libmesh_assert(!_sub_domain_ids.size());
    libmesh_assert(!_elem_subset.size());
    
    _elem_subset              = elems;
    _if_evaluate_on_all_elems = false;
}
    



void
MAST::OutputAssemblyElemOperations::
set_participating_boundaries(const std::set<libMesh::boundary_id_type>& bids) {
    
    libmesh_assert(!_bids.size());
    _bids = bids;
}
    


const std::set<const libMesh::Elem*>&
MAST::OutputAssemblyElemOperations::get_participating_elements() const {
    
    return _elem_subset;
}
    


const std::set<libMesh::subdomain_id_type>&
MAST::OutputAssemblyElemOperations::get_participating_subdomains() {
    
    return _sub_domain_ids;
}
    


const std::set<libMesh::boundary_id_type>&
MAST::OutputAssemblyElemOperations::get_participating_boundaries() {
    
    return _bids;
}
    

bool
MAST::OutputAssemblyElemOperations::
if_evaluate_for_element(const libMesh::Elem& elem) const {

    if (_if_evaluate_on_all_elems)
        return true;
    else if (_elem_subset.count(&elem))
        return true;
    else if (_sub_domain_ids.count(elem.subdomain_id()))
        return true;
    else
        return false;
}
    
    
    

