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

#ifndef __mast__second_order_newmark_transient_solver__
#define __mast__second_order_newmark_transient_solver__

// MAST includes
#include "solver/transient_solver_base.h"


namespace MAST {
    
    
    /*!
     *    This class implements the Newmark solver for solution of a
     *    second-order ODE.
     *
     */
    class SecondOrderNewmarkTransientSolver:
    public MAST::TransientSolverBase {
    public:
        SecondOrderNewmarkTransientSolver();
        
        virtual ~SecondOrderNewmarkTransientSolver();
        
        /*!
         *    \f$ \beta \f$ parameter used by this solver.
         */
        Real beta;

        /*!
         *    \f$ \gamma \f$ parameter used by this solver.
         */
        Real gamma;

        /*!
         *    update the transient velocity based on the current solution
         */
        virtual void update_velocity(libMesh::NumericVector<Real>& vec,
                                     const libMesh::NumericVector<Real>& sol);
        
        /*!
         *    update the transient acceleration based on the current solution
         */
        virtual void update_acceleration(libMesh::NumericVector<Real>& vec,
                                         const libMesh::NumericVector<Real>& sol);
        
        /*!
         *    update the transient sensitivity velocity based on the
         *    current sensitivity solution
         */
        virtual void update_sensitivity_velocity(libMesh::NumericVector<Real>& vel,
                                                 const libMesh::NumericVector<Real>& sol);
        
        /*!
         *    update the transient sensitivity acceleration based on the
         *    current sensitivity solution
         */
        virtual void update_sensitivity_acceleration(libMesh::NumericVector<Real>& acc,
                                                     const libMesh::NumericVector<Real>& sol);
        

        /*!
         *    update the perturbation in transient velocity based on the
         *    current perturbed solution
         */
        virtual void
        update_delta_velocity(libMesh::NumericVector<Real>& vel,
                              const libMesh::NumericVector<Real>& sol);
        
        /*!
         *    update the perturbation in transient acceleration based on the
         *    current perturbed solution
         */
        virtual void
        update_delta_acceleration(libMesh::NumericVector<Real>& acc,
                                  const libMesh::NumericVector<Real>& sol);

    
        /*!
         *    provides the element with the transient data for calculations
         */
        virtual void
        set_element_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                         const std::vector<libMesh::NumericVector<Real>*>& sols);
        
        /*!
         *    provides the element with the sensitivity of transient data for
         *    calculations
         */
        virtual void
        extract_element_sensitivity_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                                         const std::vector<libMesh::NumericVector<Real>*>& sols,
                                         std::vector<RealVectorX>& local_sols);

        /*!
         *    provides the element with the transient data for calculations
         */
        virtual void
        set_element_perturbed_data
        (const std::vector<libMesh::dof_id_type>& dof_indices,
         const std::vector<libMesh::NumericVector<Real>*>& sols);
        
        
        /*!
         *   performs the element calculations over \p elem, and returns
         *   the element vector and matrix quantities in \p mat and
         *   \p vec, respectively. \p if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void
        elem_calculations(bool if_jac,
                          RealVectorX& vec,
                          RealMatrixX& mat);
        
        /*!
         *   performs the element calculations over \p elem, and returns
         *   the element vector quantity in \p vec. The vector quantity only
         *   include the \f$ [J] \{dX\} f$ components, so the inherited classes
         *   must ensure that no component of constant forces (traction/body
         *   forces/etc.) are added to this vector.
         */
        virtual void
        elem_linearized_jacobian_solution_product(RealVectorX& vec);
        
        /*!
         *   performs the element sensitivity calculations over \p elem,
         *   and returns the element residual sensitivity in \p vec .
         */
        virtual void
        elem_sensitivity_calculations(const MAST::FunctionBase& f,
                                      RealVectorX& vec);

        
        /*!
         *   computes the contribution for this element from previous
         *   time step
         */
        virtual void
        elem_sensitivity_contribution_previous_timestep(const std::vector<RealVectorX>& prev_sols,
                                                        RealVectorX& vec);
        
        /*!
         *   performs the element shape sensitivity calculations over \p elem,
         *   and returns the element residual sensitivity in \p vec .
         */
        virtual void
        elem_shape_sensitivity_calculations(const MAST::FunctionBase& f,
                                            RealVectorX& vec);
        
        /*!
         *   performs the element topology sensitivity calculations over \p elem,
         *   and returns the element residual sensitivity in \p vec .
         */
        virtual void
        elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                               const MAST::LevelSetIntersection& intersect,
                                               const MAST::FieldFunction<RealVectorX>& vel,
                                               RealVectorX& vec);

        /*!
         *   calculates \f$ d ([J] \{\Delta X\})/ dX  \f$ over \p elem,
         *   and returns the matrix in \p vec .
         */
        virtual void
        elem_second_derivative_dot_solution_assembly(RealMatrixX& mat) {
            libmesh_assert(false); // to be implemented
        }

    protected:
        
    };
    
}

#endif // __mast__second_order_newmark_transient_solver__
