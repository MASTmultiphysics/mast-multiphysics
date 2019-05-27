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

#ifndef __mast_level_set_interface_dof_handler_h__
#define __mast_level_set_interface_dof_handler_h__


// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/numeric_vector.h"


namespace MAST {
    
    // Forward declerations
    class SystemInitialization;
    template <typename ValType> class FieldFunction;
    class LevelSetIntersection;
    class NonlinearImplicitAssemblyElemOperations;
    
    class LevelSetInterfaceDofHandler {
        
    public:
        
        LevelSetInterfaceDofHandler();
        
        virtual ~LevelSetInterfaceDofHandler();
        
        void init(const MAST::SystemInitialization& sys_init,
                  MAST::LevelSetIntersection& intersection,
                  MAST::FieldFunction<Real>& phi);

        MAST::FieldFunction<Real>&
        get_level_set_function() {
            libmesh_assert(_phi);
            return *_phi;
        }
        
        /*!
         *  @returns true if the element is on the interface and some nodes
         *  of this element should be factored out of the system instead of
         *  being added to the system residual vector and Jacobian matrix.
         */
        bool if_factor_element(const libMesh::Elem& elem) const;
        
        /*!
         *   identifies which rows in the element residual vector and
         *   rows/columns in the jacobian matrix correspond to the material
         *   and void. This is done under the assumption that the dofs are
         *   stored in an variable major sequence. For example, an element with
         *   four nodes and two variales will have the following sequence
         *   \f$ \{ u_1, u_2, u_3, u_4, v_1, v_2, v_3, v_4 \} \f$.
         */
        void partition_local_elem_rows(const libMesh::Elem& elem,
                                       std::vector<libMesh::dof_id_type>& material_dofs,
                                       std::vector<libMesh::dof_id_type>& void_dofs);

        
        /*!
         *   fills the \p material_dofs and \p void_dofs with the dofs_ids
         *   in the global system corresponding to these dofs in the element.
         */
        void partition_global_elem_rows(const libMesh::Elem& elem,
                                        std::vector<libMesh::dof_id_type>& material_dofs,
                                        std::vector<libMesh::dof_id_type>& void_dofs);


        /*!
         *   updates the components of the solution vector in \p elem_sol
         *   for the void domain using the stored solution for this element.
         */
        void solution_of_factored_element(const libMesh::Elem& elem,
                                          RealVectorX& elem_sol);
        

        /*!
         *   a wrapper around the second element_factored_jacobian.
         *   This takes the jacobian \p jac and factors it over the material
         *   and void domains and returns the following matrix in
         *   \p jac_factored_uu : \f$ J_{uu} - J_{uf} J_{ff}^{-1} J_{fu} \f$ .
         *   The \p *_dof_ids are the local rows of the element vector/matrix
         *   quantities that are in the material domain.
         */
        void element_factored_jacobian(const libMesh::Elem& elem,
                                       const RealMatrixX& jac,
                                       std::vector<libMesh::dof_id_type>& material_dof_ids,
                                       RealMatrixX& jac_factored_uu);

        /*!
         *    factorizes the residual and jacobian into the components for
         *    the dofs on material nodes. The factored Jacobian is defined as
         *    \p jac_factored_uu : \f$ J_{uu} - J_{uf} J_{ff}^{-1} J_{fu} \f$.
         *    The factored residual is defined as \p res_factored_u :
         *    \f$  f_{u} - J_{uf} J_{ff}^{-1} f_{f} \f$. The \p material_dof_ids
         *    are the rows in the element quantity that belong to the material
         *    domain.
         */
        void
        element_factored_residual_and_jacobian(const libMesh::Elem& elem,
                                               const RealMatrixX& jac,
                                               const RealVectorX& res,
                                               std::vector<libMesh::dof_id_type>& material_dof_ids,
                                               RealMatrixX& jac_factored_uu,
                                               RealVectorX& res_factored_u);

        /*!
         *   This takes the jacobian \p jac and factors it over the material
         *   and void domains and returns the following matrix in
         *   \p jac_factored_uu : \f$ J_{uu} - J_{uf} J_{ff}^{-1} J_{fu} \f$ .
         *   The \p *_dof_ids are the local rows of the element vector/matrix
         *   quantities that are in the material/void domain.
         */
        void element_factored_jacobian(const libMesh::Elem& elem,
                                       const RealMatrixX& jac,
                                       std::vector<libMesh::dof_id_type>& material_dof_ids,
                                       std::vector<libMesh::dof_id_type>& void_dof_ids,
                                       RealMatrixX& jac_uu,
                                       RealMatrixX& jac_uf,
                                       RealMatrixX& jac_fu,
                                       RealMatrixX& jac_ff,
                                       RealMatrixX& jac_factored_uu);

        void update_factored_element_solution(const libMesh::Elem& elem,
                                              const RealMatrixX& res,
                                              const RealMatrixX& jac,
                                              const RealMatrixX& sol,
                                              const RealMatrixX& dsol,
                                              RealVectorX& updated_sol);

        
    protected:

        void _compute_matrix_inverse(const RealMatrixX& mat,
                                     RealMatrixX& mat_inv);
        
        const MAST::SystemInitialization*           _sys_init;
        MAST::FieldFunction<Real>*                  _phi;
        
        std::map<const libMesh::Elem*, RealVectorX> _elem_sol;
        
        /*!
         *   map of nodes on each element that will be added as independent
         *   dofs.
         */
        std::map<const libMesh::Elem*, std::set<const libMesh::Node*>> _elem_void_nodes;

        /*!
         *   map of elements that each node will provide a new dof for
         */
        std::map<const libMesh::Node*, std::set<const libMesh::Elem*>> _void_node_elems;
        
        /*!
         *   new dof ids for each elem/old_dof pair.
         */
        std::map<const libMesh::Elem*, std::map<libMesh::dof_id_type, libMesh::dof_id_type>> _dof_ids;
    };
}


#endif // __mast_level_set_interface_dof_handler_h__
