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
#include "level_set/interface_dof_handler.h"
#include "level_set/level_set_intersection.h"
#include "level_set/material_patch.h"
#include "base/system_initialization.h"
#include "base/field_function_base.h"
#include "base/nonlinear_system.h"

// libMesh includes
#include "libmesh/dof_map.h"


MAST::LevelSetInterfaceDofHandler::LevelSetInterfaceDofHandler():
_sys_init   (nullptr),
_phi        (nullptr) {
    
}



MAST::LevelSetInterfaceDofHandler::~LevelSetInterfaceDofHandler() {
    
}


void
MAST::LevelSetInterfaceDofHandler::init(const MAST::SystemInitialization& sys_init,
                                        MAST::LevelSetIntersection& intersection,
                                        MAST::FieldFunction<Real>& phi) {
  
    libmesh_assert(!_sys_init);
    
    _sys_init = &sys_init;
    _phi      = &phi;
    
    const MAST::NonlinearSystem &system  = sys_init.system();
    const libMesh::DofMap       &dof_map = system.get_dof_map();
    const libMesh::MeshBase     &mesh    = system.get_mesh();

    std::set<const libMesh::Node*> negative_phi_nodes;
    std::vector<const libMesh::Elem*>
    elems_to_factor,
    negative_phi_elems;
    elems_to_factor.reserve(mesh.n_local_elem()),
    negative_phi_elems.reserve(mesh.n_local_elem());

    // iterate over all the elements identify if elements on the interface
    // need to factor their stiffness matrices
    libMesh::MeshBase::const_element_iterator
    el     =  mesh.active_local_elements_begin(),
    end_el =  mesh.active_local_elements_end();

    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        intersection.init(phi, *elem, system.time);
        if (intersection.if_elem_has_boundary()) {
            
            intersection.get_nodes_on_negative_phi(negative_phi_nodes);
            
            std::set<const libMesh::Node*>::const_iterator
            nd_it   = negative_phi_nodes.begin(),
            nd_end  = negative_phi_nodes.end();
            
            for ( ; nd_it != nd_end; nd_it++) {
                
                MAST::MaterialPatch patch;
                patch.init(*elem, **nd_it, phi, system.time);
                std::set<const libMesh::Elem*>
                patch_elems_to_factor = patch.get_elems_to_factor();
                std::set<const libMesh::Elem*>::const_iterator
                el_it  = patch_elems_to_factor.begin(),
                el_end = patch_elems_to_factor.end();
                
                for ( ; el_it != el_end; el_it++)
                    elems_to_factor.push_back(*el_it);
            }
        }
        
        negative_phi_nodes.clear();
        intersection.clear();
    }
    
    // now create a unique set from the vector of elements to facotr
    std::set<const libMesh::Elem*> elems_to_factor_set(elems_to_factor.begin(),
                                                       elems_to_factor.end());
    
    elems_to_factor.clear();
    
    // now, create an entry in the map for each one of these elements
    std::set<const libMesh::Elem*>::const_iterator
    el_it   = elems_to_factor_set.begin(),
    el_end  = elems_to_factor_set.end();
    
    std::vector<libMesh::dof_id_type> dof_ids;
    std::pair<const libMesh::Elem*, RealVectorX> elem_sol_pair;
    
    for ( ; el_it != el_end; el_it++) {
        
        dof_map.dof_indices(*el_it, dof_ids);
        elem_sol_pair.first = *el_it;
        elem_sol_pair.second = RealVectorX::Zero(dof_ids.size());
        _elem_sol.insert(elem_sol_pair);
    }
}


bool
MAST::LevelSetInterfaceDofHandler::if_factor_element(const libMesh::Elem& elem) const {
    
    libmesh_assert(_sys_init);
    
    return _elem_sol.count(&elem);
}


void
MAST::LevelSetInterfaceDofHandler::
partition_local_elem_rows(const libMesh::Elem& elem,
                          std::vector<libMesh::dof_id_type>& material_dofs,
                          std::vector<libMesh::dof_id_type>& void_dofs) {
    
    
    const MAST::NonlinearSystem &system  = _sys_init->system();
    Real val = 0.;
    unsigned int
    nvars   = _sys_init->n_vars(),
    n_nodes = elem.n_nodes();
    
    // currently only implemented for quad4 and Lagrange functions
    libmesh_assert_equal_to(elem.type(), libMesh::QUAD4);
    // the number of dofs is assumed to be a multiple of the number of variables.
    // currently, it is assumed that the Lagrange shape functions are used,
    // so that each node has dofs for each variable.
    libmesh_assert_equal_to(system.variable_type(0).family, libMesh::LAGRANGE);

    // the dofs are assumed to be sequenced in variable major format
    // So, we loop over the nodes and if a node is in the void domain, all
    // dofs on that will be factored out.
    material_dofs.reserve(nvars*elem.n_nodes());
    void_dofs.reserve(nvars*elem.n_nodes());

    for (unsigned int i=0; i<elem.n_nodes(); i++) {
        
        (*_phi)(*elem.node_ptr(i), system.time, val);
        if (val < 0.) {
            
            // mark all dofs for this node to be in void
            for (unsigned int j=0; j<nvars; j++)
                void_dofs.push_back(j*n_nodes+i);
        }
        else {
            
            // mark all dofs for this node to be in material
            for (unsigned int j=0; j<nvars; j++)
                material_dofs.push_back(j*n_nodes+i);
        }
    }
    
    // now sort the dofs
    std::sort(void_dofs.begin(), void_dofs.end());
    std::sort(material_dofs.begin(), material_dofs.end());
}


void
MAST::LevelSetInterfaceDofHandler::
partition_global_elem_rows(const libMesh::Elem& elem,
                           std::vector<libMesh::dof_id_type>& material_dofs,
                           std::vector<libMesh::dof_id_type>& void_dofs) {
    
    // get the local dof partitioning
    this->partition_local_elem_rows(elem, material_dofs, void_dofs);
    
    //
    // now replace the dof ids with the global ids
    //
    const MAST::NonlinearSystem &system  = _sys_init->system();
    const libMesh::DofMap       &dof_map = system.get_dof_map();
    std::vector<libMesh::dof_id_type>
    dof_ids;
    
    dof_map.dof_indices(&elem, dof_ids);
    
    for (unsigned int i=0; i<material_dofs.size(); i++)
        material_dofs[i] = dof_ids[material_dofs[i]];
    
    for (unsigned int i=0; i<void_dofs.size(); i++)
        void_dofs[i]     = dof_ids[void_dofs[i]];
}


void
MAST::LevelSetInterfaceDofHandler::solution_of_factored_element(const libMesh::Elem& elem,
                                                                RealVectorX& elem_sol) {
    
    // make sure that the element has been identified for factorization
    std::map<const libMesh::Elem*, RealVectorX>::const_iterator
    it   = _elem_sol.find(&elem),
    end  = _elem_sol.end();
    
    libmesh_assert(it != end);
    
    std::vector<libMesh::dof_id_type>
    material_dofs,
    void_dofs;

    const RealVectorX
    &free_sol = it->second;
    
    this->partition_local_elem_rows(elem, material_dofs, void_dofs);
    
    // overwirte the void dofs from the stored data
    for (unsigned int i=0; i<void_dofs.size(); i++)
        elem_sol(void_dofs[i]) = free_sol(void_dofs[i]);
}



void
MAST::LevelSetInterfaceDofHandler::
element_factored_jacobian(const libMesh::Elem& elem,
                          const RealMatrixX& jac,
                          std::vector<libMesh::dof_id_type>& material_dof_ids,
                          RealMatrixX& jac_factored_uu) {
    
    std::vector<libMesh::dof_id_type>
    void_dof_ids;
    
    RealMatrixX
    jac_uu,
    jac_uf,
    jac_fu,
    jac_ff;
    
    this->element_factored_jacobian(elem,
                                    jac,
                                    material_dof_ids,
                                    void_dof_ids,
                                    jac_uu,
                                    jac_uf,
                                    jac_fu,
                                    jac_ff,
                                    jac_factored_uu);
    
}



void
MAST::LevelSetInterfaceDofHandler::
element_factored_residual_and_jacobian(const libMesh::Elem& elem,
                                       const RealMatrixX& jac,
                                       const RealVectorX& res,
                                       std::vector<libMesh::dof_id_type>& material_dof_ids,
                                       RealMatrixX& jac_factored_uu,
                                       RealVectorX& res_factored_u) {
    
    std::vector<libMesh::dof_id_type>
    void_dof_ids;
    
    RealMatrixX
    jac_uu,
    jac_uf,
    jac_fu,
    jac_ff,
    jac_ff_inv;
    
    this->element_factored_jacobian(elem,
                                    jac,
                                    material_dof_ids,
                                    void_dof_ids,
                                    jac_uu,
                                    jac_uf,
                                    jac_fu,
                                    jac_ff,
                                    jac_factored_uu);

    unsigned int
    nu   = material_dof_ids.size(),
    nf   = void_dof_ids.size();

    RealVectorX
    res_f           = RealVectorX::Zero(nf);
    
    res_factored_u  = RealVectorX::Zero(nu);
    
    for (unsigned int i=0; i<material_dof_ids.size(); i++)
        res_factored_u(i)  = res(material_dof_ids[i]);
    
    for (unsigned int i=0; i<void_dof_ids.size(); i++)
        res_f(i)           = res(void_dof_ids[i]);
    
    _compute_matrix_inverse(jac_ff, jac_ff_inv);
    
    res_factored_u -= jac_uf * jac_ff_inv * res_f;
}



void
MAST::LevelSetInterfaceDofHandler::
element_factored_jacobian(const libMesh::Elem& elem,
                          const RealMatrixX& jac,
                          std::vector<libMesh::dof_id_type>& material_dof_ids,
                          std::vector<libMesh::dof_id_type>& void_dof_ids,
                          RealMatrixX& jac_uu,
                          RealMatrixX& jac_uf,
                          RealMatrixX& jac_fu,
                          RealMatrixX& jac_ff,
                          RealMatrixX& jac_factored_uu) {
    

    material_dof_ids.clear();
    void_dof_ids.clear();
    
    this->partition_local_elem_rows(elem, material_dof_ids, void_dof_ids);
    
    unsigned int
    nu = material_dof_ids.size(),
    nf = void_dof_ids.size();
    
    jac_uu = RealMatrixX::Zero(nu, nu);
    jac_uf = RealMatrixX::Zero(nu, nf);
    jac_fu = RealMatrixX::Zero(nf, nu);
    jac_ff = RealMatrixX::Zero(nf, nf);
    
    RealMatrixX
    jac_ff_inv;
    
    //
    // partition the matrices
    //
    for (unsigned int i=0; i<nu; i++) {
        for (unsigned int j=0; j<nu; j++)
            jac_uu(i,j) = jac(material_dof_ids[i], material_dof_ids[j]);

        for (unsigned int j=0; j<nf; j++)
            jac_uf(i,j) = jac(material_dof_ids[i],     void_dof_ids[j]);
    }

    for (unsigned int i=0; i<nf; i++) {
        for (unsigned int j=0; j<nu; j++)
            jac_fu(i,j) = jac(void_dof_ids[i], material_dof_ids[j]);
        
        for (unsigned int j=0; j<nf; j++)
            jac_ff(i,j) = jac(void_dof_ids[i],     void_dof_ids[j]);
    }
    
    _compute_matrix_inverse(jac_ff, jac_ff_inv);
    
    jac_factored_uu = jac_uu - jac_uf * jac_ff_inv * jac_fu;
}


void
MAST::LevelSetInterfaceDofHandler::
update_factored_element_solution(const libMesh::Elem& elem,
                                 const RealMatrixX& res,
                                 const RealMatrixX& jac,
                                 const RealMatrixX& sol,
                                 const RealMatrixX& dsol,
                                 RealVectorX& updated_sol) {
    
    std::map<const libMesh::Elem*, RealVectorX>::iterator
    it  = _elem_sol.find(&elem),
    end = _elem_sol.end();
    
    libmesh_assert(it != end);
    
    //
    //   [Juu  Juf] {dxu} = - {f_u}
    //   [Jfu  Jff] {dxf} = - {f_f}
    //   so,
    //   {dxf} = inv(Jff) (-{f_f} - Jfu {dxu})
    //

    std::vector<libMesh::dof_id_type>
    material_dof_ids,
    void_dof_ids;
    
    this->partition_local_elem_rows(elem, material_dof_ids, void_dof_ids);
    
    // we should not be capturing elements that are purely in void or material
    libmesh_assert_greater(material_dof_ids.size(), 0);
    libmesh_assert_greater(    void_dof_ids.size(), 0);

    unsigned int
    nu = material_dof_ids.size(),
    nf = void_dof_ids.size();
    
    RealMatrixX
    jac_fu = RealMatrixX::Zero(nf, nu),
    jac_ff = RealMatrixX::Zero(nf, nf),
    jac_ff_inv;
    
    RealVectorX
    f_u    = RealVectorX::Zero(nu),
    f_f    = RealVectorX::Zero(nf),
    dx_u   = RealVectorX::Zero(nu),
    dx_f   = RealVectorX::Zero(nf);
    

    //
    // partition the matrices
    //
    for (unsigned int i=0; i<nf; i++) {
        
        for (unsigned int j=0; j<nu; j++)
            jac_fu(i,j) = jac(void_dof_ids[i], material_dof_ids[j]);
        
        for (unsigned int j=0; j<nf; j++)
            jac_ff(i,j) = jac(void_dof_ids[i],     void_dof_ids[j]);
    
        f_f(i) = res(void_dof_ids[i]);
    }

    for (unsigned int i=0; i<nu; i++)
        dx_u(i)= dsol(material_dof_ids[i]);

    _compute_matrix_inverse(jac_ff, jac_ff_inv);
    
    dx_f = jac_ff_inv * (-f_f - jac_fu * dx_u);
    
    updated_sol = sol;
    for (unsigned int i=0; i<nf; i++)
        updated_sol(void_dof_ids[i]) += dx_f(i);
    

    //
    // now update the solution in the map
    //
    it->second = updated_sol;
}


void
MAST::LevelSetInterfaceDofHandler::_compute_matrix_inverse(const RealMatrixX& mat,
                                                           RealMatrixX& mat_inv) {
    
    const unsigned int
    m = mat.rows(),
    n = mat.cols();
    
    libmesh_assert_equal_to(m, n);
    
    Eigen::FullPivLU<RealMatrixX> solver(mat);
    RealVectorX
    rhs = RealVectorX::Zero(m);
    
    mat_inv = RealMatrixX::Zero(m, m);
    
    for (unsigned int i=0; i<m; i++) {
        
        rhs  = RealVectorX::Zero(m);
        rhs(i) = 1.;
        
        mat_inv.col(i) = solver.solve(rhs);
    }
}

